#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>

double generate_sig2 (double *Y_Xbeta, double *Z, double *gamma, double *Rinv, int n, int p, int s, 
		double phi2, double a, double b)
{
/* Note that A = Y_Xbeta - Z %*% gamma */
  double alpha=-1.0, beta1=1.0, working;
  double a1, b1, norm_gamma, norm_A, rsig2, *tmp_prod;
  int inc=1;
  char notrans='n';

  tmp_prod = (double *) calloc(n, sizeof(double));
  if (tmp_prod==NULL) {
    Rprintf("Malloc failed for tmp_prod.\n");
    exit(EXIT_FAILURE);
  }

  a1 = (n+s)/2 + a;

/* NOTE THAT AFTER THIS Y_XBETA HAS BEEN OVERWRITTEN 
   AND IS NOW EQUAL TO (Y_XBETA - Z %*% GAMMA) !!!!!!!*/
  F77_CALL(dgemv)(&notrans, &n, &s, &alpha, Z, &n, gamma, &inc, &beta1, Y_Xbeta, &inc);
  alpha=1.0;
  beta1=0.0;
  F77_CALL(dgemv)(&notrans, &n, &n, &alpha, Rinv, &n, Y_Xbeta, &inc, &beta1, tmp_prod, &inc);
  norm_A = F77_CALL(ddot)(&n, Y_Xbeta, &inc, tmp_prod, &inc);
  free(tmp_prod);
/* t(Y - X %*% beta) %*% Rinv %*% (Y - X %*% beta)  COMPUTED*/

  working = F77_CALL(dnrm2)(&s, gamma, &inc);
 
  norm_gamma = R_pow( working, 2.0);
  b1 = (norm_A + norm_gamma/phi2)/2.0 + b;
  GetRNGstate();
  rsig2 = 1.0/rgamma(a1, 1.0/b1);  
  PutRNGstate();
  return rsig2;
}
