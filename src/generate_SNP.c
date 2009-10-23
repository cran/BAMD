#include<R.h>
#include<Rmath.h>
#include<math.h>
#include<R_ext/BLAS.h>
#include"SNP.h"

void generate_SNP (double *beta, double *gamma, double *X, double *Y, double *Z, 
		int n, int p, int s, double sig2, int missing_index, double **exact_missing,
 		int *row_missing, int *col_missing, int num_missing, double *Zprob, 
		double *R, int SNPprior)
{
  double *Z1, tmp_Y_Xbeta, num[3], p1[3], prior[3], working;
  int inc=1, i, three=3, inc2=n;
  char notrans='n', trans='t';

  i=missing_index;
  Z1 = (double *) malloc(sizeof(double) * s);
  if (Z1==NULL) {
    Rprintf("Malloc failed for Z1\n");
    exit(EXIT_FAILURE);
  }

  F77_CALL(dcopy)(&s, Z+row_missing[i], &n, Z1, &inc);
  tmp_Y_Xbeta = *(Y + row_missing[i]) - F77_CALL(ddot)(&p, X + row_missing[i], &n, beta, &inc);

  /*If informative priors used, FBLAS */
  if(SNPprior==1) {
    F77_CALL(dcopy)(&three, Zprob + i, &num_missing, prior, &inc);

      *(Z1+col_missing[i])=ALL1;
      num[0]= prior[0] * exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));

      *(Z1+col_missing[i])=ALL2;
      num[1]= prior[1] * exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));

      *(Z1+col_missing[i])=ALL3;
      num[2]= prior[2] * exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));
  }

  /*If noninformative priors used, FBLAS */
  else {
    F77_CALL(dcopy)(&s, Z+row_missing[i], &n, Z1, &inc);
    tmp_Y_Xbeta = *(Y + row_missing[i]) - F77_CALL(ddot)(&p, X + row_missing[i], &n, beta, &inc);

      *(Z1+col_missing[i])=ALL1;
      num[0]= exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));

      *(Z1+col_missing[i])=ALL2;
      num[1]= exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));

      *(Z1+col_missing[i])=ALL3;
      num[2]= exp(R_pow(tmp_Y_Xbeta - F77_CALL(ddot)(&s, Z1, &inc, gamma, &inc), 2.0)/ (-2 * sig2));
  }

   working=num[0]+num[1]+num[2];
   p1[0]=num[0]/working;
   p1[1]=num[1]/working;

  GetRNGstate();
  working = unif_rand();
  PutRNGstate();

   if(working<=p1[0])
	**(exact_missing+i) = ALL1;
   else if (working<=(p1[0]+p1[1]))
	**(exact_missing+i) = ALL2;
   else 
	**(exact_missing+i) = ALL3;
  
  free(Z1);
}
