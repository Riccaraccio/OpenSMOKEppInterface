#include <stdio.h>
#include <math.h>
#include "OpenSMOKE_Interface.h"

#define kinfolder "../kinetics/skeletal/methanol/kinetics"

int main () {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);

  const int NS = OpenSMOKE_NumberOfSpecies();
  const char* key_species[] = {"CH3OH", "O2"};

  // CH3OH/air with a small radical pool, at 1500 K and 1 atm
  double x[NS], c[NS], R[NS], Rred[NS];
  for (int i=0; i<NS; i++)
    x[i] = 0.;
  x[OpenSMOKE_IndexOfSpecies ("CH3OH")] = 0.12;
  x[OpenSMOKE_IndexOfSpecies ("O2")]    = 0.18;
  x[OpenSMOKE_IndexOfSpecies ("N2")]    = 0.68;
  x[OpenSMOKE_IndexOfSpecies ("H")]     = 0.005;
  x[OpenSMOKE_IndexOfSpecies ("OH")]    = 0.005;
  x[OpenSMOKE_IndexOfSpecies ("HO2")]   = 0.005;
  x[OpenSMOKE_IndexOfSpecies ("H2O")]   = 0.005;

  double T = 1500., P = 101325.;
  for (int i=0; i<NS; i++)
    c[i] = P/(8314.462618*T)*x[i];

  OpenSMOKE_GasProp_SetTemperature (T);
  OpenSMOKE_GasProp_SetPressure (P);
  OpenSMOKE_GasProp_ReactionRates (c);
  OpenSMOKE_GasProp_FormationRates (R);
  double QR = OpenSMOKE_GasProp_HeatRelease (R);

  // Number of important species and reactions for increasing eps.
  // With eps = 0, the reduced rates must be equal to the full rates.
  double eps_list[] = {0., 1e-2, 5e-2, 1e-1, 3e-1};
  for (double eps : eps_list) {
    double Tmax = 1e5;
    OpenSMOKE_DRG_Init (key_species, 2, &eps, &Tmax, 1, 300., 0.);
    OpenSMOKE_DRG_Analysis (T, P, c);
    OpenSMOKE_DRG_ReducedFormationRates (c, Rred);

    double dR = 0.;
    for (int i=0; i<NS; i++)
      dR = fmax (dR, fabs (Rred[i] - R[i]));

    fprintf (stderr, "eps = %g: ns = %d, nr = %d\n", eps,
             OpenSMOKE_DRG_NumberImportantSpecies(),
             OpenSMOKE_DRG_NumberImportantReactions());
    if (eps == 0.)
      fprintf (stderr, "eps = 0: max |dR| = %g, dQR = %g\n",
               dR, OpenSMOKE_GasProp_HeatRelease (Rred) - QR);
  }

  // Below the threshold temperature: only the key species, no reactions
  OpenSMOKE_DRG_Analysis (250., P, c);
  fprintf (stderr, "T = 250 K: ns = %d, nr = %d\n",
           OpenSMOKE_DRG_NumberImportantSpecies(),
           OpenSMOKE_DRG_NumberImportantReactions());

  OpenSMOKE_Clean();
  return 0;
}
