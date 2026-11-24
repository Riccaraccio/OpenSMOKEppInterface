// Note: this file is included inside namespace OpenSMOKE by OpenSMOKE_Jacobian.h,
//       so standard headers must be included there, not here.

GasJacobian::GasJacobian(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), kineticsMap_(kineticsMap)
{
	NS_ = thermodynamicsMap_.NumberOfSpecies();
	NR_ = kineticsMap_.NumberOfReactions();

	is_supported_ = (kineticsMap_.NumberOfCABRReactions() == 0);

	ones_.resize(NS_, 1.);
	c_fd_.resize(NS_, 0.);
	rf_.resize(NR_, 0.);
	rb_.resize(NR_, 0.);
	uKeq_.resize(NR_, 0.);

	// Reaction orders: the matrices are NR x NS (column-major),
	// so it.row() is the reaction and it.col() is the species
	reactant_species_.resize(NR_);
	reactant_orders_.resize(NR_);
	product_species_.resize(NR_);
	product_orders_.resize(NR_);
	{
		const Eigen::SparseMatrix<double>& orders = kineticsMap_.stoichiometry().reactionorders_matrix_reactants();
		for (int k = 0; k < orders.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(orders, k); it; ++it)
			{
				reactant_species_[it.row()].push_back(it.col());
				reactant_orders_[it.row()].push_back(it.value());
			}
	}
	{
		// Only the reversible reactions have product orders
		const Eigen::SparseMatrix<double>& orders = kineticsMap_.stoichiometry().reactionorders_matrix_products();
		for (int k = 0; k < orders.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(orders, k); it; ++it)
			{
				product_species_[it.row()].push_back(it.col());
				product_orders_[it.row()].push_back(it.value());
			}
	}

	// Net stoichiometric coefficients nu = products - reactants (dense NR x NS, then lists)
	{
		Eigen::MatrixXd nu(NR_, NS_);
		nu.setZero();

		const Eigen::SparseMatrix<double>& products = kineticsMap_.stoichiometry().stoichiometric_matrix_products();
		for (int k = 0; k < products.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(products, k); it; ++it)
				nu(it.row(), it.col()) += it.value();

		const Eigen::SparseMatrix<double>& reactants = kineticsMap_.stoichiometry().stoichiometric_matrix_reactants();
		for (int k = 0; k < reactants.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(reactants, k); it; ++it)
				nu(it.row(), it.col()) -= it.value();

		nu_species_.resize(NR_);
		nu_.resize(NR_);
		for (unsigned int r = 0; r < NR_; r++)
			for (unsigned int i = 0; i < NS_; i++)
				if (nu(r, i) != 0.)
				{
					nu_species_[r].push_back(i);
					nu_[r].push_back(nu(r, i));
				}
	}
}

double GasJacobian::Power(const double c, const double a)
{
	if (a == 1.)
		return c;
	else if (a == 2.)
		return c*c;
	else if (a == 3.)
		return c*c*c;
	else if (a == 0.5)
		return std::sqrt(c);
	else
		return std::pow(c, a);
}

double GasJacobian::DerivativeOfPower(const double c, const double a)
{
	// c^(a-1) is infinite at c = 0 for a < 1 and not defined for c < 0
	// with a non-integer order: use max(c, c_floor) in this factor only
	const double c_floor = 1.e-20;

	if (a == 1.)
		return 1.;
	else if (a == 2.)
		return 2.*c;
	else if (a == 3.)
		return 3.*c*c;
	else
		return a*std::pow(std::max(c, c_floor), a - 1.);
}

void GasJacobian::AddToColumn(const unsigned int r, const unsigned int j, const double d, double* dRdc)
{
	double* column = dRdc + NS_*j;
	for (unsigned int i = 0; i < nu_species_[r].size(); i++)
		column[nu_species_[r][i]] += nu_[r][i]*d;
}

void GasJacobian::AddOrders(const unsigned int r, const double k, const double* c,
	const std::vector<unsigned int>& species, const std::vector<double>& orders, double* dRdc)
{
	// d(k prod_m c_m^a_m)/dc_j = k a_j c_j^(a_j-1) prod_{m != j} c_m^a_m
	// The product of the other species avoids a division by c_j (c_j can be 0 or negative)
	for (unsigned int n = 0; n < species.size(); n++)
	{
		double d = k*DerivativeOfPower(c[species[n]], orders[n]);
		for (unsigned int m = 0; m < species.size(); m++)
			if (m != n)
				d *= Power(c[species[m]], orders[m]);
		AddToColumn(r, species[n], d, dRdc);
	}
}

void GasJacobian::Compute(const double* c, double* dRdc)
{
	// 1. kb/kf from a call at c = 1: then rf = kf and rb = kf*kb/kf
	kineticsMap_.ReactionRates(ones_.data());
	kineticsMap_.GetForwardReactionRates(rf_.data());
	std::fill(rb_.begin(), rb_.end(), 0.);			// only the reversible reactions are written
	kineticsMap_.GetBackwardReactionRates(rb_.data());
	for (unsigned int r = 0; r < NR_; r++)
		uKeq_[r] = (rf_[r] != 0.) ? rb_[r]/rf_[r] : 0.;	// rf = 0 if the DRG turned the reaction off

	// 2. Rates at c (cTot is kept constant in the derivatives)
	double cTot = 0.;
	for (unsigned int i = 0; i < NS_; i++)
		cTot += c[i];

	kineticsMap_.ReactionRates(c, cTot);
	kineticsMap_.GetForwardReactionRates(rf_.data());
	std::fill(rb_.begin(), rb_.end(), 0.);
	kineticsMap_.GetBackwardReactionRates(rb_.data());
	const std::vector<double>& kf = kineticsMap_.KArrheniusModified();	// with third-body and fall-off factors

	std::fill(dRdc, dRdc + NS_*NS_, 0.);

	// 3. Reaction orders: d(rf - rb)/dc_j at constant kf
	for (unsigned int r = 0; r < NR_; r++)
	{
		if (kf[r] == 0.)
			continue;
		AddOrders(r, kf[r], c, reactant_species_[r], reactant_orders_[r], dRdc);
		if (uKeq_[r] != 0.)
			AddOrders(r, -kf[r]*uKeq_[r], c, product_species_[r], product_orders_[r], dRdc);
	}

	// 4. Third-body reactions: kf is proportional to M = cTot + sum_k e_k c_k
	for (unsigned int s = 0; s < kineticsMap_.NumberOfThirdBodyReactions(); s++)
	{
		const unsigned int r = kineticsMap_.IndicesOfThirdbodyReactions()[s] - 1;
		const double M = kineticsMap_.M()[s];
		if (M == 0.)
			continue;
		const double f = (rf_[r] - rb_[r])/M;
		for (unsigned int k = 0; k < kineticsMap_.IndicesOfThirdbodySpecies()[s].size(); k++)
			AddToColumn(r, kineticsMap_.IndicesOfThirdbodySpecies()[s][k] - 1,
				f*kineticsMap_.IndicesOfThirdbodyEfficiencies()[s][k], dRdc);
	}

	// 5. Fall-off reactions: kf is proportional to the correction G(c),
	//    central finite difference of G in c_j (the efficiencies are private)
	std::copy(c, c + NS_, c_fd_.begin());
	for (unsigned int l = 0; l < kineticsMap_.NumberOfFallOffReactions(); l++)
	{
		const unsigned int r = kineticsMap_.IndicesOfFalloffReactions()[l] - 1;
		const double G = kineticsMap_.FallOffReactionsCorrection(l + 1, cTot, c);
		if (G == 0.)
			continue;
		const double f = (rf_[r] - rb_[r])/G;

		const unsigned int single = kineticsMap_.FallOffIndexOfSingleThirdbodySpecies()[l];
		if (single != 0)
		{
			// M = c_s + 1e-16
			const unsigned int j = single - 1;
			const double h = 1.e-6*(std::fabs(c[j]) + 1.e-16);
			c_fd_[j] = c[j] + h;
			const double Gp = kineticsMap_.FallOffReactionsCorrection(l + 1, cTot, c_fd_.data());
			c_fd_[j] = c[j] - h;
			const double Gm = kineticsMap_.FallOffReactionsCorrection(l + 1, cTot, c_fd_.data());
			c_fd_[j] = c[j];
			AddToColumn(r, j, f*(Gp - Gm)/(2.*h), dRdc);
		}
		else
		{
			// M = cTot + sum_k e_k c_k
			const double h = 1.e-6*cTot;
			for (unsigned int k = 0; k < kineticsMap_.FallOffIndicesOfThirdbodySpecies()[l].size(); k++)
			{
				const unsigned int j = kineticsMap_.FallOffIndicesOfThirdbodySpecies()[l][k] - 1;
				c_fd_[j] = c[j] + h;
				const double Gp = kineticsMap_.FallOffReactionsCorrection(l + 1, cTot, c_fd_.data());
				c_fd_[j] = c[j] - h;
				const double Gm = kineticsMap_.FallOffReactionsCorrection(l + 1, cTot, c_fd_.data());
				c_fd_[j] = c[j];
				AddToColumn(r, j, f*(Gp - Gm)/(2.*h), dRdc);
			}
		}
	}
}
