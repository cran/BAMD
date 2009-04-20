#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

void generate_mvn(double *mean, double *covm, double *sim, int n) {
	int info, inc=1, i;
	char uplo='L';
	char notrans='n', nonunit='n';
	double alpha=1.0;
	
	F77_CALL(dpotrf)(&uplo, &n, covm, &n, &info);	

  	GetRNGstate();
	for(i=0; i<n; i++)
		sim[i] = rnorm(0.0, 1.0);
  	PutRNGstate();

	F77_CALL(dtrmv)(&uplo, &notrans, &nonunit, &n, covm, &n, sim, &inc);
	F77_CALL(daxpy)(&n, &alpha, mean, &inc, sim, &inc);
}
