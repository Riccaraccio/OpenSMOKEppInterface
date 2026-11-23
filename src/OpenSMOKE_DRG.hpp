// Copied from LaminarSMOKEpp by A. Cuoci
// Note: this file is included inside namespace OpenSMOKE by OpenSMOKE_DRG.h,
//       so standard headers must be included there, not here.

DRG::DRG(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
    OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, const std::vector<std::string>& key_species, 
    const std::vector<double>& temperatures, const double thresholdTemperature,
    const double thresholdFormationRate, const std::vector<double>& epsilons):
  thermodynamicsMap_(thermodynamicsMap), kineticsMap_(kineticsMap)
{
  NS_ = thermodynamicsMap_.NumberOfSpecies();
  NR_ = kineticsMap_.NumberOfReactions();

  rNet_.resize(NR_); rNet_.setZero();
  rReduced_.resize(NR_); rReduced_.setZero();
  r_.resize(NS_, NS_);
  numerator_.resize(NS_);
  denominator_.resize(NS_);
  Rp_.resize(NS_);
  Rd_.resize(NS_);

  number_important_species_ = NS_;
  number_unimportant_reactions_ = 0;
  number_important_reactions_ = NR_;

  threshold_temperature_ = 300.;
  formation_rate_threshold_ = 1e-12;

  // Build a full matrix of net stoichiometric coefficients nu = nuB - nuF
  Eigen::MatrixXd nu_(NR_, NS_);						
  {
    // Be careful: eigen vectors and matrices are 0-index based
    // Be careful: if the kinetic scheme is large, this matrix, since it is full, can be very memory expensive
    //             Example: 10^3 species, 10^4 reactions = size of the matrix 10^7 elements!
    //             This is the reason why we store stoichiometric matrices in sparse format.
    //             Of course te advantage of having a full matrix, is that you access the elements directly, without
    //             using iterators and pointers, as reported above
    nu_.setZero();

    // Loop over all the reactions (product side)
    for (int k = 0; k < kineticsMap_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
    {
      // Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
      for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap_.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
      {
        nu_(it.row(), it.col()) += it.value();
      }
    }

    // Loop over all the reactions (product side)
    for (int k = 0; k < kineticsMap_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
    {
      // Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
      for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMap_.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
      {
        nu_(it.row(), it.col()) -= it.value();
      }
    }
  }

  // Build the delta matrix (dense matrix) used by the DRG method
  Eigen::MatrixXi delta_(NR_, NS_);					
  {
    for (unsigned int i = 0; i < NR_; i++)
    {
      for (unsigned int j = 0; j < NS_; j++)
        delta_(i, j) = (nu_(i, j) == 0.) ? 0 : 1;
    }
  }

  // Sparse Algebra
  {			
    // Sparse net stoichiometric matrix
    {
      nu_sparse_.resize(NS_,NR_);

      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.reserve(NR_*4);
      for (unsigned int i = 0; i < NR_; i++)
        for (unsigned int j = 0; j < NS_; j++)
        {
          if ( nu_(i,j) != 0.)
            tripletList.push_back(T(j,i,nu_(i,j)));
        }

      nu_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    // Sparse delta matrix, 1 means the species is involved in the reaction, (NR x NS)
    delta_sparse_.resize(NS_,NR_);
    {
      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.reserve(NR_*4);
      for (unsigned int i = 0; i < NR_; i++)
        for (unsigned int j = 0; j < NS_; j++)
        {
          if ( delta_(i,j) != 0 )
            tripletList.push_back(T(j,i,1.));
        }

      delta_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    // Nu times delta
    {
      nu_times_delta_.resize(NS_);

      for (unsigned int k = 0; k < NS_; k++)
      {
        nu_times_delta_[k].resize(NS_, NR_);

        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(NR_*4);
        for (unsigned int i = 0; i < NR_; i++)
          for (unsigned int j = 0; j < NS_; j++)
          {
            const double prod = nu_(i,k) * delta_(i,j);
            if ( prod != 0.)
              tripletList.push_back(T(j,i, prod));
          }

        nu_times_delta_[k].setFromTriplets(tripletList.begin(), tripletList.end());
      }
    }
  }

  // Initialize DRG object
  {
    // Read list of ke-species
    Eigen::VectorXi drgListSpecies;
    drgListSpecies.resize(key_species.size());
    for (int i=0;i<drgListSpecies.size();i++)
      drgListSpecies(i) = thermodynamicsMap_.IndexOfSpecies(key_species[i])-1;

    // Set key species
    SetKeySpecies(drgListSpecies);

    // Read list of temperature/epsilon pairs
    epsilons_ = epsilons;
    temperatures_ = temperatures;

    // Read threshold temperature
    threshold_temperature_ = thresholdTemperature; 

    // Read threshold formation rate
    formation_rate_threshold_ = thresholdFormationRate;
  }
}

DRG::~DRG() {}

void DRG::SetKeySpecies(const std::vector<std::string>& names_key_species)
{
  index_key_species_.resize(names_key_species.size());
  for (unsigned int i = 0; i < names_key_species.size(); i++)
    index_key_species_(i) = thermodynamicsMap_.IndexOfSpecies(names_key_species[i]) - 1;
}

void DRG::SetKeySpecies(const Eigen::VectorXi& key_species)
{
  index_key_species_ = key_species;
}

void DRG::Analysis(const double T, const double P_Pa, const double* c)
{
  // Look for threshold (epsilon)
  double epsilon = 1e-2;
  for(unsigned int i=0;i<temperatures_.size();i++)
    if(T < temperatures_[i]) {
      epsilon = epsilons_[i];
      break;
    }

  // Perform DRG
  PairWiseErrorMatrix(T, P_Pa, c);
  ParsePairWiseErrorMatrix(epsilon, T);
}

void DRG::PairWiseErrorMatrix(const double T, const double P_Pa, const double* c)
{
  // Now we know T, P and composition. 
  // We have to pass those data to the thermodynamic and kinetic maps
  thermodynamicsMap_.SetTemperature(T);
  thermodynamicsMap_.SetPressure(P_Pa);
  kineticsMap_.SetTemperature(T);
  kineticsMap_.SetPressure(P_Pa);

  // Now we can calculate (internally) the reaction rates concentrations are needed
  kineticsMap_.ReactionRates(c);
  kineticsMap_.GiveMeReactionRates(rNet_.data());				// [kmol/m3/s]

  // Now we calculate the production and consumption rates of species
  kineticsMap_.ProductionAndDestructionRates(Rp_.data(), Rd_.data());	// [kmol/m3/s]

  // Calculate the pair-wise error matrix
  r_.setZero();

  if (T >= threshold_temperature_)
  {
    // Denominator			
    denominator_.setConstant(0.);		
    for (int k=0; k<nu_sparse_.outerSize(); ++k)
    {
      const double rnet = rNet_(k);
      for (Eigen::SparseMatrix<double>::InnerIterator it(nu_sparse_,k); it; ++it)
      {
        denominator_(it.row()) += std::fabs(it.value() * rnet);
      }
    }

    // Numerator
    for (unsigned int i = 0; i < NS_; i++)
    {
      // If the production and consumption rates of a species are below a threshold,
      // we can put the pair-wise error matrix equal to zero, because the reaction is negligible
      if ( Rp_(i) <= formation_rate_threshold_ && Rd_(i) <= formation_rate_threshold_ )
      {
        // nothing to do
      }
      else
      {
        numerator_.setConstant(0.);
        for (int k=0; k<nu_times_delta_[i].outerSize(); ++k)
        {
          const double rnet = rNet_(k);
          for (Eigen::SparseMatrix<double>::InnerIterator it(nu_times_delta_[i],k); it; ++it)
          {
            numerator_(it.row()) += std::fabs(it.value() * rnet);
          }
        }

        for (unsigned int j = 0; j < NS_; j++)
          r_(i,j) = numerator_(j)/(1.e-64+denominator_(i));
      }
    }	
  }
}

void DRG::ParsePairWiseErrorMatrix(const double epsilon, const double T)
{
  // Perform the analysis only in case of sufficiently high temperature
  if (T >= threshold_temperature_)
  {
    std::vector<bool> vector_important_species(NS_);
    std::vector<bool> vector_important_reactions(NR_);

    // Reset important important reactions
    vector_important_species.assign(NS_,false);
    vector_important_reactions.assign(NR_,true);

    // Initialize the queue with key-species
    std::queue <int> Q;
    for ( int i=0; i<index_key_species_.size(); i++)
    {
      Q.push(index_key_species_(i));
      vector_important_species[index_key_species_(i)] = true;
    }

    // DFS with rAB
    while (!Q.empty())
    {
      for (unsigned int k=0; k<NS_; k++)
      {
        if (vector_important_species[k] == false)                    
        {
          if (r_(Q.front(), k) > epsilon)
          {                            
            vector_important_species[k] = true;
            Q.push(k);
          }
        }
      }
      Q.pop();
    }	

    // Important reactions
    for (int k=0; k<delta_sparse_.outerSize(); ++k)
    {
      for (Eigen::SparseMatrix<double>::InnerIterator it(delta_sparse_,k); it; ++it)
      {
        if (vector_important_species[it.row()] == false)
        {
          vector_important_reactions[k] = false;
        }
      }
    }

    // Count important species and reactions
    number_important_species_ = std::count (vector_important_species.begin(), vector_important_species.end(), true);
    number_unimportant_reactions_ = std::count (vector_important_reactions.begin(), vector_important_reactions.end(), false);
    number_important_reactions_ = NR_ - number_unimportant_reactions_;

    // Vector containing the indices of important species (zero-based)
    {
      indices_important_species_.resize(number_important_species_);
      unsigned int count = 0;
      for(unsigned int k = 0; k < NS_; k++)
        if (vector_important_species[k] == true)
        {
          indices_important_species_(count) = k;
          count++;
        }
    }

    // Vector containing the indices of unimportant reactions (zero-based)
    {
      indices_unimportant_reactions_.resize(number_unimportant_reactions_);
      indices_important_reactions_.resize(number_important_reactions_);

      unsigned int count_unimportant = 0;
      unsigned int count_important = 0;
      for(unsigned int k = 0; k < NR_; k++)
        if (vector_important_reactions[k] == false)
        {
          indices_unimportant_reactions_(count_unimportant) = k;
          count_unimportant++;
        }
        else
        {
          indices_important_reactions_(count_important) = k;
          count_important++;
        }
    }

    // Vector containing the indices of important reactions steps (zero-based)
    {
      important_reactions_steps_.resize(number_important_reactions_);

      if (number_important_reactions_ > 0)
      {
        important_reactions_steps_(0) = indices_important_reactions_(0);
        for(unsigned int k = 1; k < number_important_reactions_; k++)
          important_reactions_steps_(k) = indices_important_reactions_(k) - indices_important_reactions_(k-1);
      }
    }
  }
  else
  {
    // Count important species and reactions
    number_important_species_ = index_key_species_.size();
    number_unimportant_reactions_ = NR_;
    number_important_reactions_ = 0;

    // Vector containing the indices of important species (zero-based)
    {
      indices_important_species_ = index_key_species_;
    }

    // Vector containing the indices of unimportant reactions (zero-based)
    {
      indices_unimportant_reactions_.resize(number_unimportant_reactions_);
      indices_important_reactions_.resize(number_important_reactions_);

      for(unsigned int k = 0; k < NR_; k++)
        indices_unimportant_reactions_(k) = k;

      important_reactions_steps_.resize(number_important_reactions_);
    }
  }
}

void DRG::ReducedFormationRates(const double* c, double* R)
{
  kineticsMap_.ReactionRates(c);                 // all constants, normal cache
  kineticsMap_.GiveMeReactionRates(rReduced_.data());
  for (int i = 0; i < indices_unimportant_reactions_.size(); i++)
    rReduced_(indices_unimportant_reactions_(i)) = 0.;
  kineticsMap_.stoichiometry().FormationRatesFromReactionRates(R, rReduced_.data());
}

