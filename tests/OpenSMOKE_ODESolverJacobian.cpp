// Test of OpenSMOKE_ODESolver_Jacobian(): the Robertson problem with its
// analytic Jacobian, against OpenSMOKE_ODESolver() (finite-difference
// Jacobian). The two solutions must agree to the tolerance of the solver
// (rtol 1.2e-5, atol 1e-10), the solver must call the Jacobian, and it
// must make fewer calls of the equations. The log (stderr) has "ok"/"fail"
// only; the values go to stdout.

#include <stdio.h>
#include <math.h>
#include "OpenSMOKE_Interface.h"

typedef struct {
  double k1, k2, k3;
  int nrhs, njac;
} UserData;

static void robertson (const double * y, const double t, double * dy, void * args) {
  UserData * data = (UserData *)args;
  data->nrhs++;
  dy[0] = -data->k1*y[0] + data->k2*y[1]*y[2];
  dy[1] =  data->k1*y[0] - data->k2*y[1]*y[2] - data->k3*y[1]*y[1];
  dy[2] =  data->k3*y[1]*y[1];
}

// J[i + 3*j] = df_i/dy_j (column-major)
static void robertson_jacobian (const double * y, const double t, double * J, void * args) {
  UserData * data = (UserData *)args;
  data->njac++;
  const double k1 = data->k1, k2 = data->k2, k3 = data->k3;
  J[0 + 3*0] = -k1;  J[0 + 3*1] =  k2*y[2];               J[0 + 3*2] =  k2*y[1];
  J[1 + 3*0] =  k1;  J[1 + 3*1] = -k2*y[2] - 2.*k3*y[1];  J[1 + 3*2] = -k2*y[1];
  J[2 + 3*0] =  0.;  J[2 + 3*1] =  2.*k3*y[1];            J[2 + 3*2] =  0.;
}

int main (void) {
  OpenSMOKE_InitODESolver ();

  const int neq = 3;
  const double dts[] = {1., 40., 4e5};
  for (double dt : dts) {
    UserData fd = {0.04, 1e4, 3e7, 0, 0}, an = fd;
    double yfd[] = {1., 0., 0.}, yan[] = {1., 0., 0.};

    OpenSMOKE_ODESolver (&robertson, neq, dt, yfd, &fd);
    OpenSMOKE_ODESolver_Jacobian (&robertson, &robertson_jacobian, neq, dt, yan, &an);

    int agree = 1;
    for (int i = 0; i < neq; i++)
      if (fabs (yan[i] - yfd[i]) > 10.*(1.2e-5*fabs (yfd[i]) + 1e-10))
        agree = 0;

    fprintf (stderr, "dt = %g: solutions agree: %s, Jacobian called: %s, "
             "fewer calls of the equations: %s\n", dt, agree ? "ok" : "fail",
             an.njac > 0 && fd.njac == 0 ? "ok" : "fail",
             an.nrhs < fd.nrhs ? "ok" : "fail");
    printf ("dt = %g: FD y = %.8e %.8e %.8e, calls %d\n", dt, yfd[0], yfd[1],
            yfd[2], fd.nrhs);
    printf ("dt = %g: AN y = %.8e %.8e %.8e, calls %d, Jacobians %d\n", dt,
            yan[0], yan[1], yan[2], an.nrhs, an.njac);
  }

  OpenSMOKE_CleanODESolver ();
  return 0;
}
