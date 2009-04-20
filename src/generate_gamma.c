#include<R.h>
#include "gibbs.h"
#include<Rmath.h>
#include<R_ext/BLAS.h>
#include<R_ext/Lapack.h>

void generate_gamma (double *beta, double *X, double *Y, double *Y_Xbeta,
		double *Z, double *rgamma, double *Rinv, int n, int p, int s, 
		double phi2, double sig2)
{
  double *mean, *covm, alpha=1.0, beta1=0.0;
  double *A, *Ainv, *B, *tmp_Zt_Rinv;
  int *ipiv, tmp1;
  int inc=1, tmp;
  char notrans='n', trans='t';

  
  mean = (double *) calloc(s, sizeof(double));
  covm = (double *) calloc(s * s, sizeof(double));
  A = (double *) calloc(s * s, sizeof(double));
  B = (double *) calloc(s, sizeof(double));
  Ainv = (double *) calloc(s * s, sizeof(double));
  tmp_Zt_Rinv = (double *) calloc(s * n, sizeof(double));
  ipiv = (int *) malloc(s * sizeof(int));
  
  if ((mean==NULL)||(covm==NULL)) {
    Rprintf("Malloc failed for mean, covm.\n");
    exit(EXIT_FAILURE);
  }
  if ((B==NULL)||(A==NULL)||(ipiv==NULL)||(Ainv==NULL)) {
    Rprintf("Malloc failed for Ainv, B, ipiv, A.\n");
    exit(EXIT_FAILURE);
  }

/* construct A which will eventually be		*/
/*  A = Zt %*% R_inv %*% Z + I/phi2  	*/
  for(tmp=0;tmp<s;tmp++) *(A + tmp*s +tmp)=1.0/phi2;
  for(tmp=0;tmp<s;tmp++) *(Ainv + tmp*s +tmp)=1.0;

  F77_CALL(dgemm)(&trans, &notrans, &s, &n , &n, &alpha, Z, &n, Rinv, &n, &beta1, tmp_Zt_Rinv, &s);
  beta1=1.0;
  F77_CALL(dgemm)(&notrans, &notrans, &s, &s , &n, &alpha, tmp_Zt_Rinv, &s, Z, &n, &beta1, A, &s);
/* A constructed */

/* TODO: Compute the inverse here */
/* compute Ainv */
  F77_CALL(dgesv)(&s, &s, A, &s, ipiv, Ainv, &s, &tmp1);
  if (tmp1!=0) {
    Rprintf("Matrix inversion unsuccessful (Zt %*% R_inv %*% Z + I/phi2)!\n");
    exit(EXIT_FAILURE);
  }
  free(A);
/* Ainv computed */
  
/* compute (Y - X %*% beta) */
  F77_CALL(dcopy)(&n, Y, &inc, Y_Xbeta, &inc);

  alpha=-1.0;

  F77_CALL(dgemv)(&notrans, &n, &p, &alpha, X, &n, beta, &inc, &beta1, Y_Xbeta, &inc);
/* (Y - X %*%  beta) computed */ 

/* compute B= Zt %*% Rinv %*% (Y - X %*% beta) */
  alpha=1.0;
  beta1=0.0;

  F77_CALL(dgemv)(&notrans, &s, &n, &alpha, tmp_Zt_Rinv, &s, Y_Xbeta, &inc, &beta1, B, &inc);
  free(tmp_Zt_Rinv);
/* B = Zt %*% Rinv %*% (Y - X %*% beta) computed */


/* compute mean vector of gamma */ 
  F77_CALL(dgemv)(&notrans, &s, &s, &alpha, Ainv, &s, B, &inc, &beta1, mean, &inc);
  free(B);
/* mean vector of gamma computed */ 

/* compute covariance matrix of gamma */ 
/* which is given by (sig2 * Ainv) */
  tmp = s*s;
  F77_CALL(dcopy)(&tmp, Ainv, &inc, covm, &inc);
  F77_CALL(dscal)(&tmp, &sig2, covm, &inc);
/* covariance matrix of gamma computed */ 

  generate_mvn(mean, covm, rgamma, s);

  free(mean);
  free(covm);
  free(Ainv);
}
