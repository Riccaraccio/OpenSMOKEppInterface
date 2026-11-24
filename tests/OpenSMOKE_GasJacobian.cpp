// Test of OpenSMOKE_GasProp_FormationRatesJacobian() against a central
// finite difference of the formation rates.
//
// The Jacobian holds cTot constant. OpenSMOKE_GasProp_ReactionRates() sets
// cTot = sum of c, so the test makes the difference along v = e_j - e_k
// (sum of v = 0, cTot constant). This gives J(:,j) - J(:,k). The test does
// not see an error that is the same in all the columns: at constant P, the
// chain rule to dY removes such a term (sum_k dc_k/dY_j = 0).
//
// The test uses the C API only. If it includes the OpenSMOKE++ headers, the
// globals of those headers exist two times, and the program aborts at exit.
//
// The log (stderr) has "ok"/"fail" only, because the errors change in the
// last digits between machines. The errors and the times go to stdout.

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "OpenSMOKE_Interface.h"

// Formation rates at c + h*(e_j - e_k), with cTot = sum of c
static void Rates (const double* c, int n, int j, int k, double h, double* R) {
  std::vector<double> cp (c, c + n);
  cp[j] += h;
  cp[k] -= h;
  OpenSMOKE_GasProp_ReactionRates (cp.data());
  OpenSMOKE_GasProp_FormationRates (R);
}

// Central difference along e_j - e_k, with a Richardson extrapolation:
// D = (4 D(h/2) - D(h))/3. The round-off error of a difference is
// eps*(gross rate)/h: with h = 1e-6*c_j the columns of the inert species
// have errors of 1e-3. The extrapolation removes the h^2 error, so a large
// step (1e-2) is possible.
static void FiniteDifference (const double* c, int n, int j, int k,
                              double h, double* D) {
  std::vector<double> Rp (n), Rm (n), Rp2 (n), Rm2 (n);
  Rates (c, n, j, k, h, Rp.data());
  Rates (c, n, j, k, -h, Rm.data());
  Rates (c, n, j, k, h/2., Rp2.data());
  Rates (c, n, j, k, -h/2., Rm2.data());
  for (int i = 0; i < n; i++)
    D[i] = (4.*(Rp2[i] - Rm2[i])/h - (Rp[i] - Rm[i])/(2.*h))/3.;
}

// Largest difference in a column, relative to the largest entry of that column
static double ColumnError (const double* J, const double* c, int n) {
  int k = 0;
  double cTot = 0.;
  for (int i = 0; i < n; i++) {
    cTot += c[i];
    if (c[i] > c[k])
      k = i;
  }

  double err = 0.;
  std::vector<double> D (n);
  for (int j = 0; j < n; j++) {
    if (j == k)
      continue;
    const double h = 1e-2*fmax (fabs (c[j]), 1e-3*cTot);
    FiniteDifference (c, n, j, k, h, D.data());
    double norm = 0., diff = 0.;
    for (int i = 0; i < n; i++) {
      const double Jv = J[i + n*j] - J[i + n*k];
      norm = fmax (norm, fabs (D[i]));
      diff = fmax (diff, fabs (Jv - D[i]));
    }
    if (norm > 0.)
      err = fmax (err, diff/norm);
    else if (diff > 0.)
      err = HUGE_VAL;
  }
  return err;
}

static void TestMechanism (const char* kinfolder, double T, double P) {
  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder);
  const int n = OpenSMOKE_NumberOfSpecies();
  fprintf (stderr, "%s: %d species\n", kinfolder, n);

  OpenSMOKE_GasProp_SetTemperature (T);
  OpenSMOKE_GasProp_SetPressure (P);

  // Three states: all species equal; all species different; some species
  // at 0 and one below 0 (the RHS does not clip the mass fractions)
  std::vector<double> x (n), c (n), J (n*n);
  for (int state = 0; state < 3; state++) {
    for (int i = 0; i < n; i++)
      x[i] = state == 0 ? 1. : 1.5 + sin (1. + i);
    if (state == 2) {
      for (int i = 0; i < n; i += 3)
        x[i] = 0.;
      x[n - 2] = -1e-10;
    }
    double xTot = 0.;
    for (int i = 0; i < n; i++)
      xTot += x[i];
    for (int i = 0; i < n; i++)
      c[i] = P/(8314.462618*T)*x[i]/xTot;

    int status = OpenSMOKE_GasProp_FormationRatesJacobian (c.data(), J.data());
    double err = ColumnError (J.data(), c.data(), n);

    const double tol = 1e-5;
    fprintf (stderr, "  state %d: status %d, column error < %g: %s\n", state,
             status, tol, err < tol ? "ok" : "fail");
    printf ("%s state %d: column error %.2e\n", kinfolder, state, err);
  }

  // Time for each call (same T as the call before, as in the RHS)
  const int ncalls = 200;
  struct timespec t0, t1;
  clock_gettime (CLOCK_MONOTONIC, &t0);
  for (int k = 0; k < ncalls; k++)
    OpenSMOKE_GasProp_FormationRatesJacobian (c.data(), J.data());
  clock_gettime (CLOCK_MONOTONIC, &t1);
  printf ("%s: %.1f us for each Jacobian\n", kinfolder,
          1e6*((t1.tv_sec - t0.tv_sec) + 1e-9*(t1.tv_nsec - t0.tv_nsec))/ncalls);

  OpenSMOKE_Clean();
}

int main () {
  TestMechanism ("../kinetics/skeletal/hydrogen/kinetics", 1500., 101325.);
  TestMechanism ("../kinetics/skeletal/methanol/kinetics", 1500., 101325.);
  TestMechanism ("../kinetics/biomass/Solid-gas-88/kinetics", 1500., 101325.);
  return 0;
}
