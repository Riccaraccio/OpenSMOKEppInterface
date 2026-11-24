// Copied from LaminarSMOKEpp by A. Cuoci
#pragma once

#include "maps/Maps_CHEMKIN"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <cmath>

namespace OpenSMOKE {
class DRG {
public:

	/**
	*@brief Constructor
	*@param thermodynamicsMap thermodynamic map
	*@param kineticsMap kinetics map
	*@param key_species names of the key (target) species, always retained as important
	*@param temperatures upper temperature bounds [K] of the intervals for selecting epsilon (ascending order)
	*@param thresholdTemperature temperature [K] below which the DRG analysis is skipped
	*       (only key species are kept and all reactions are flagged as unimportant)
	*@param thresholdFormationRate production/destruction rate threshold [kmol/m3/s] below which
	*       a species is considered inactive (its row of the pair-wise error matrix is set to zero)
	*@param epsilons DRG thresholds, one per entry of temperatures: epsilons[i] is used for the first i
	*       such that T < temperatures[i]; if T is above all bounds, epsilon = 1e-2 is used
	*/
  DRG(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
    OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap, const std::vector<std::string>& key_species, 
    const std::vector<double>& temperatures, const double thresholdTemperature,
    const double thresholdFormationRate, const std::vector<double>& epsilons);

	/**
	*@brief Default destructor
	*/
	virtual ~DRG();

	/**
	*@brief Sets the key-species (i.e. the target species which are considered important)
	*@param names_key_species vector containing the names of target or key species (0-index based)
	*/
	void SetKeySpecies(const std::vector<std::string>& names_key_species);

	/**
	*@brief Sets the key-species (i.e. the target species which are considered important)
	*@param names_key_species vector containing the indices of target or key species (0-index based)
	*/
	void SetKeySpecies(const Eigen::VectorXi& indices_key_species);

	/**
	*@brief Performs the DRG analysis for the given conditions
	*@param T temperature in K
	*@param P_Pa pressure in Pa
	*@param c vector of concentrations in kmol/m3
	*/
	void Analysis(const double T, const double P_Pa, const double* c);

	/**
	*@brief Returns the number of key species
	*/
	unsigned int number_key_species() const { return index_key_species_.size(); }

	/**
	*@brief Returns the number of important species
	*/
	unsigned int number_important_species() const { return number_important_species_; }

	/**
	*@brief Returns the number of important reactions
	*/
	unsigned int number_important_reactions() const { return number_important_reactions_; }

	/**
	*@brief Returns the indices of important species (zero-based)
	*/
	const Eigen::VectorXi& indices_important_species() const { return indices_important_species_; }

	/**
	*@brief Returns the indices of important reactions (zero-based)
	*/
	const Eigen::VectorXi& indices_important_reactions() const { return indices_important_reactions_; }

	/**
	*@brief Returns the indices of unimportant reactions (zero-based)
	*/
	const Eigen::VectorXi& indices_unimportant_reactions() const { return indices_unimportant_reactions_; }

	/**
	*@brief Returns the indices of important reactions steps (zero-based)
	*/
	const Eigen::VectorXi& important_reactions_steps() const { return important_reactions_steps_; }

  /**
   *@brief Formation rates [kmol/m3/s] of the reduced mechanism, from the last Analysis
   *       Set T and P on the maps before you call this function.
   *       The function computes the kinetic constants of the important reactions only
   *       (as laminarSMOKE++). The kinetic map then holds 1e-64 for the constants of the
   *       unimportant reactions, and its cache marks the constants as valid. Call
   *       RestoreKineticConstants() before a full use of the kinetic map.
   */
  void ReducedFormationRates(const double* c, double* R);

  /**
   *@brief Makes the kinetic map compute all the kinetic constants again at its next use,
   *       if ReducedFormationRates() changed them. The map keeps the temperature T.
   *@param T temperature [K] that the kinetic map must hold after the call
   */
  void RestoreKineticConstants(const double T);

private:

	/**
	*@brief Builds the pair wise error matrix
	*@param T temperature in K
	*@param P_Pa pressure in Pa
	*@param c vector of concentrations in kmol/m3 
	*@param epsilon threshold
	*/
	void PairWiseErrorMatrix(const double T, const double P_Pa, const double* c);

  /**
   *@brief Analyzes the pair wise error matrix and calculates (internally) the important_species vector
   */
  void ParsePairWiseErrorMatrix(const double epsilon, const double T);

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;		/**< reference to the thermodynamic map */
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;				/**< reference to the kinetics map */

	Eigen::VectorXi index_key_species_;					/**< indices of target or key species (zero-based) */
  Eigen::VectorXd rReduced_;

	Eigen::VectorXd rNet_;							/**< vector containing the net reaction rates */
	unsigned int NR_;							/**< total number of reactions */
	unsigned int NS_;							/**< total number of species */

	unsigned int number_important_species_;					/**< current number of important species */
	unsigned int number_important_reactions_;				/**< number of unimportant reactions */
	unsigned int number_unimportant_reactions_;				/**< number of unimportant reactions */

	Eigen::VectorXi indices_important_reactions_;				/**< indices of important reactions (zero-based) */
	Eigen::VectorXi indices_unimportant_reactions_;				/**< indices of unimportant reactions (zero-based) */
	Eigen::VectorXi indices_important_species_;				/**< indices of important species (zero-based) */
	Eigen::VectorXi important_reactions_steps_;				/**< indices of important reactions steps (zero-based) */

	Eigen::MatrixXd					r_;			/**< full matrix, pair wise error matrix, (NS x NS) */
	Eigen::SparseMatrix<double> 			nu_sparse_;		/**< sparse matrix containing the net stoichiometric coefficients, (NS x NR) */
	std::vector< Eigen::SparseMatrix<double> >	nu_times_delta_;	/**< vector of sparse matrices (NS x NS) */
	Eigen::SparseMatrix<double> 			delta_sparse_;		/**< delta matrix (NS x NR) */

  bool reduced_constants_;   /**< true if the kinetic map holds the constants of the reduced mechanism */

  std::vector<double>  epsilons_;
  std::vector<double>  temperatures_;
	double threshold_temperature_;
	double formation_rate_threshold_;

	Eigen::VectorXd Rp_;
	Eigen::VectorXd Rd_;
	Eigen::VectorXd numerator_;
	Eigen::VectorXd denominator_;
};

#include "OpenSMOKE_DRG.hpp"
}
