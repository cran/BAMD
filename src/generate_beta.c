#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include"gibbs.h"

void generate_beta (double *rbeta, double *XtX_inv_Xt, double *XtXinv, double *Y, 
		double *Z, double *gamma, int n, int p, int s, double sig2)
{
  double *mean, *covm, *tmp_Y_Zgam, alpha=-1.0, beta1=1.0;
  int inc=1, tmp;
  char notrans='n', trans='t';

  mean = (double *) calloc(p, sizeof(double));
  covm = (double *) calloc(p * p, sizeof(double));

  if ((mean==NULL)||(covm==NULL)) {
    Rprintf("Malloc failed for mean, covm.\n");
    exit(EXIT_FAILURE);
  }

/* compute (Y - Z %*% gamma)				*/
  tmp_Y_Zgam = (double *) calloc(n, sizeof(double));
  if (tmp_Y_Zgam==NULL) {
    Rprintf("Malloc failed for tmp_Y_Zgam\n");
    exit(EXIT_FAILURE);
  }

  F77_CALL(dcopy)(&n, Y, &inc, tmp_Y_Zgam, &inc);
  F77_CALL(dgemv)(&notrans, &n, &s, &alpha, Z, &n, gamma, &inc, &beta1, tmp_Y_Zgam, &inc);
/* (Y - Z %*%  gamma) computed */

/* compute mean vector of beta, which is given by (XtX_inv %*% Xt %*% (Y - Z %*% gamma)) */
  beta1=0.0;
  alpha=1.0;

  F77_CALL(dgemv)(&notrans, &p, &n, &alpha, XtX_inv_Xt, &p, tmp_Y_Zgam, &inc, &beta1, mean, &inc);
/* mean vector of beta computed */ 

/* compute covariance matrix of beta, which is given by (sig2 * XtXinv) */
  tmp = p*p;
  F77_CALL(dcopy)(&tmp, XtXinv, &inc, covm, &inc);
  F77_CALL(dscal)(&tmp, &sig2, covm, &inc);
/* covariance matrix of beta vector computed */

  generate_mvn(mean, covm, rbeta, p);
  
  free(tmp_Y_Zgam);
  free(mean);
  free(covm);
}
