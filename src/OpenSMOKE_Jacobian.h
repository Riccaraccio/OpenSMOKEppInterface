#pragma once

#include "maps/Maps_CHEMKIN"
#include <vector>
#include <algorithm>
#include <cmath>

namespace OpenSMOKE {

/**
*@brief Analytical Jacobian of the formation rates of the gas species with respect to the
*       concentrations, at constant T, P and total concentration cTot.
*       dR_i/dc_j = sum_r nu_ri * D(r,j), where D(r,j) = d(rf_r - rb_r)/dc_j.
*       The class uses only the public data of KineticsMap_CHEMKIN:
*       - the effective forward constant kf_r is KArrheniusModified() after ReactionRates(c)
*       - the ratio kb_r/kf_r comes from a call of ReactionRates() with c = 1 for all species
*       - the fall-off efficiencies are private: the derivatives of the fall-off correction
*         use a central finite difference, one for each (reaction, species) pair
*       Supported: elementary, reversible (thermodynamic and explicit), third body,
*       fall-off (Lindemann, Troe, SRI), PLOG, Chebyshev, non-integer orders (FORD).
*       Not supported: CABR reactions (is_supported() is false). The extended fall-off and
*       extended PLOG reactions are not supported, and the class cannot detect them,
*       because their counts are private.
*/
class GasJacobian
{
public:

	/**
	*@brief Constructor: makes the lists of the reaction orders and of the net
	*       stoichiometric coefficients
	*@param thermodynamicsMap thermodynamic map
	*@param kineticsMap kinetics map
	*/
	GasJacobian(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Computes the Jacobian. Set T and P on the maps before you call this function.
	*       The function calls ReactionRates() of the kinetic map two times: the map holds
	*       the rates at c after the call.
	*@param c concentrations [kmol/m3]
	*@param dRdc dense Jacobian, column-major, NS x NS: dRdc[i+NS*j] = dR_i/dc_j [1/s]
	*/
	void Compute(const double* c, double* dRdc);

	/**
	*@brief Returns false if the mechanism has a reaction type that the class does not support
	*/
	bool is_supported() const { return is_supported_; }

private:

	/**
	*@brief c^a, as StoichiometricMap::ProductOfConcentrations (regularized for a < 1)
	*/
	static double Power(const double c, const double a);

	/**
	*@brief d(c^a)/dc (regularized for a < 1, finite at c = 0)
	*/
	static double DerivativeOfPower(const double c, const double a);

	/**
	*@brief Adds D(r,j) = d to the column j of the Jacobian: dRdc(i,j) += nu_ri * d
	*/
	void AddToColumn(const unsigned int r, const unsigned int j, const double d, double* dRdc);

	/**
	*@brief Adds the derivatives of one side of reaction r (k * prod c^order)
	*/
	void AddOrders(const unsigned int r, const double k, const double* c,
		const std::vector<unsigned int>& species, const std::vector<double>& orders, double* dRdc);

private:

	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_;	/**< reference to the thermodynamic map */
	OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_;			/**< reference to the kinetics map */

	unsigned int NR_;						/**< total number of reactions */
	unsigned int NS_;						/**< total number of species */
	bool is_supported_;						/**< false if the mechanism has CABR reactions */

	std::vector< std::vector<unsigned int> > reactant_species_;	/**< species of the reactant orders of each reaction (zero-based) */
	std::vector< std::vector<double> >       reactant_orders_;	/**< reactant orders of each reaction */
	std::vector< std::vector<unsigned int> > product_species_;	/**< species of the product orders of each reaction (zero-based) */
	std::vector< std::vector<double> >       product_orders_;	/**< product orders of each reaction */
	std::vector< std::vector<unsigned int> > nu_species_;		/**< species with a net stoichiometric coefficient (zero-based) */
	std::vector< std::vector<double> >       nu_;			/**< net stoichiometric coefficients (products - reactants) */

	std::vector<double> ones_;		/**< c = 1 for all species */
	std::vector<double> c_fd_;		/**< concentrations for the fall-off finite differences */
	std::vector<double> rf_;		/**< forward reaction rates [kmol/m3/s] */
	std::vector<double> rb_;		/**< backward reaction rates [kmol/m3/s] */
	std::vector<double> uKeq_;		/**< kb/kf of each reaction (0 for irreversible reactions) */
};

#include "OpenSMOKE_Jacobian.hpp"
}
