#include<R.h>
#include<R_ext/BLAS.h>
#include<R_ext/Lapack.h>
#include"gibbs.h"

void gibbs_sampler(double *Y, double *X, double *Z, double *R, double *Zprob, int *missing_index, 
		int *num_missing, int *nsim, int *n, int *p, int *s, double *a, double *b, double *c, 
		double *d, double *beta_s, double *gamma_s, double *sig2_s, double *phi2_s, int *paramet, 
		int *keep, int *SNPprior)
{  
/* Variables that will be kept and re-used again and again */
  double *Xt_Rinv_X_inv, *Xt_Rinv_X_inv_Xt_Rinv, *Y_Xbeta, *Rinv, **exact_missing;
  int *row_missing, *col_missing;
  FILE *fmissing;
  

/* Disposable variables */
  double *tmp_XtRX, *tmp_XtRinv, alpha=1.0, beta=0.0;
  int i, j=0, k=0, SNPcol, loc_n, loc_p, loc_s, num_miss, parametrization;
  int loc_SNPprior;
  int *ipiv, info;
  char notrans='n', trans='t';

  loc_n = *n;
  loc_p = *p;
  loc_s = *s;
  parametrization = *paramet;

  num_miss = *num_missing;
  loc_SNPprior = *SNPprior;

/* Malloc steps (perm vars) - make sure these are free-d!*/
  Xt_Rinv_X_inv = (double *) calloc(loc_p * loc_p, sizeof(double));
  Xt_Rinv_X_inv_Xt_Rinv = (double *) calloc(loc_p * loc_n, sizeof(double));
  Y_Xbeta = (double *) calloc(loc_n, sizeof(double));
  Rinv = (double *) calloc(loc_n * loc_n, sizeof(double));
  exact_missing = malloc(num_miss * sizeof(double *));
  if ((Rinv==NULL)||(Xt_Rinv_X_inv==NULL)||(Xt_Rinv_X_inv_Xt_Rinv==NULL)||(Y_Xbeta==NULL)||(exact_missing==NULL)) {
    Rprintf("Malloc failed for XtX_inv, XtX_inv_Xt, Y_Xbeta, exact_missing\n");
    exit(EXIT_FAILURE);
  }

/* Malloc steps (temp vars) - make sure these are free-d immediately ! */
  tmp_XtRX = (double *) calloc(loc_p * loc_p, sizeof(double));
  tmp_XtRinv = (double *) calloc(loc_p * loc_n, sizeof(double));
  ipiv = (int *) calloc(loc_n, sizeof(int));
  row_missing = (int *) calloc(num_miss, sizeof(int));
  col_missing = (int *) calloc(num_miss, sizeof(int));
  if ((tmp_XtRinv==NULL)||(ipiv==NULL)||(tmp_XtRX==NULL)||(row_missing==NULL)||(col_missing==NULL)) {
    Rprintf("Malloc failed for tmp_XtX, ipiv\n");
    exit(EXIT_FAILURE);
  }

/* COMPUTE R_inverse */
  for(i=0; i<loc_n; i++) *(Rinv + i*loc_n + i)=1.0;
  F77_CALL(dgesv)(&loc_n, &loc_n, R, &loc_n, ipiv, Rinv, &loc_n, &info);
  if (info!=0) {
    Rprintf("Matrix inversion unsuccessful (R-inverse)!\n");
    exit(EXIT_FAILURE);
  }
/* Rinv COMPUTED */
  Rprintf("R-inverse computed.\n");
  

/* ASSIGN THE ARRAY OF POINTERS FIRST */
  for(i=0;i<num_miss;i++) {
    *(row_missing + i) = *(missing_index+i);
    *(col_missing + i) = *(missing_index+i+num_miss);
    if(parametrization==1) 
      *(exact_missing + i) = Z + *(row_missing + i) + (*(col_missing + i) * loc_n) ;
    else 
      *(exact_missing + i) = Z + *(row_missing + i) + (*(col_missing + i) * 2 * loc_n) ;
  }

  /* SET-UP OF (t(X) %*% Rinv %*% X)inv */
  /* First set-up the identity matrix */
  for(i=0;i<loc_p;i++) *(Xt_Rinv_X_inv + i*loc_p +i)=1.0;
  F77_CALL(dgemm)(&trans, &notrans, &loc_p, &loc_n, &loc_n, &alpha, X, &loc_n, Rinv, &loc_n, &beta, tmp_XtRinv, &loc_p);
  F77_CALL(dgemm)(&notrans, &notrans, &loc_p, &loc_p, &loc_n, &alpha, tmp_XtRinv, &loc_p, X, &loc_n, &beta, tmp_XtRX, &loc_p);

  F77_CALL(dgesv)(&loc_p, &loc_p, tmp_XtRX, &loc_p, ipiv, Xt_Rinv_X_inv, &loc_p, &info);
  if (info!=0) {
    Rprintf("Matrix inversion unsuccessful (solve(t(X) %*% Rinv %*% X))!\n");
    exit(EXIT_FAILURE);
  }
  free(tmp_XtRX);
  free(ipiv);
  Rprintf("Xt_Rinv_X_inv computed.\n");
/* (t(X) %*% Rinv %*% X)inv COMPUTED */

  /* SET-UP OF (X' %*% Rinv %*% X)inv %*% X' %*% Rinv */
  F77_CALL(dgemm)(&notrans, &notrans, &loc_p, &loc_n, &loc_p, &alpha, Xt_Rinv_X_inv, &loc_p, tmp_XtRinv, &loc_p, &beta, Xt_Rinv_X_inv_Xt_Rinv, &loc_p);
  free(tmp_XtRinv);
  /* (X' %*% Rinv %*% X)inv %*% X' %*% Rinv COMPUTED */
    
  Rprintf("Entering loop..\n");
/* THIS IS WHERE ALL THE ACTION BEGINS  */

  for(i=1;i<*nsim;i++) {
  Rprintf("generate beta, sim=%d\n", i);
  generate_beta (beta_s+i*loc_p, Xt_Rinv_X_inv_Xt_Rinv, Xt_Rinv_X_inv, Y, Z, gamma_s+(i-1)*loc_s, 
	loc_n, loc_p, loc_s, *(sig2_s+i-1)); 
  
  Rprintf("generate gamma, sim=%d\n", i); 
  generate_gamma (beta_s+i*loc_p, X, Y, Y_Xbeta, Z, gamma_s+i*loc_s, Rinv, loc_n, 
	loc_p, loc_s, *(phi2_s+i-1), *(sig2_s+i-1));
  
  Rprintf("generate sig2, sim=%d\n", i); 
  *(sig2_s+i) = generate_sig2 (Y_Xbeta, Z, gamma_s+i*loc_s, Rinv, loc_n, loc_p, loc_s, 
	*(phi2_s+i-1), *a, *b);
  
  Rprintf("generate phi2, sim=%d\n", i); 
  *(phi2_s+i) = generate_phi2 (gamma_s+i*loc_s, loc_s, *(sig2_s+i), *c, *d);
  
  Rprintf("generate SNPcol, sim=%d\n", i); 
  SNPcol= *(col_missing+j);

  if(num_miss>1) {
  while(*(col_missing+j) == SNPcol) {
/* update one missing value at a time; one SNP col at a time */
  generate_SNP (beta_s+i*loc_p, gamma_s+i*loc_s, X, Y, Z, loc_n, loc_p, loc_s, 
	*(sig2_s+i), j, exact_missing, row_missing, col_missing, num_miss, 
	Zprob, R, parametrization, loc_SNPprior);
  if(i>(*nsim - *keep)) {
			  fprintf(fmissing, "%.1f ", **(exact_missing + j));
  }
  j = (++j)%num_miss;
  }
  if(i>(*nsim - *keep)) {
	  fprintf(fmissing, "\n");
  	if(i<(*nsim-1)) 
	  fprintf(fmissing, "%d ", *(col_missing + j) + 1);
  }
  }
  else if(num_miss==1) {
  generate_SNP (beta_s+i*loc_p, gamma_s+i*loc_s, X, Y, Z, loc_n, loc_p, loc_s, 
	*(sig2_s+i), j, exact_missing, row_missing, col_missing, num_miss, 
	Zprob, R, parametrization, loc_SNPprior);
  if(i>(*nsim - *keep)) {
		  fprintf(fmissing, "%.1f ", **(exact_missing + j));
	  fprintf(fmissing, "\n");
  	if(i<(*nsim-1)) 
	  fprintf(fmissing, "%d ", *(col_missing + j) + 1);
  }
  }
  if(i==(*nsim - *keep)) { 
	  fmissing = fopen("Imputed_missing_vals", "w");
	  for(k=0;k<num_miss;k++) {
		  fprintf(fmissing, "%.1f ", **(exact_missing + k));
	  }
          if(num_miss>0) {
	  fprintf(fmissing, "\n");
	  fprintf(fmissing, "%d ", *(col_missing + j) + 1);
	  }
  }

  }
  Rprintf("Gibbs sampler loop finished!\n");
  
  fclose(fmissing);
  free(Xt_Rinv_X_inv);
  free(Xt_Rinv_X_inv_Xt_Rinv);
  free(Y_Xbeta);
  free(Rinv);
  free(exact_missing);
  free(row_missing);
  free(col_missing);
}
