#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>

double generate_phi2 (double *gamma, int s, double sig2, double c, double d)
{
  double working;
  double a2, b2, rphi2;
  int inc=1;

  a2 = s/2.0 + c;

  working =  F77_CALL(dnrm2)(&s, gamma, &inc);

  b2 = R_pow( working, 2.0 )/ (2.0*sig2) + d;
  GetRNGstate();
  rphi2 = 1.0/rgamma(a2, 1.0/b2);
  PutRNGstate();
  
  return rphi2;
}
