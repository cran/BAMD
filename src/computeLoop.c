#include<R.h>
#include<Rmath.h>
#include<math.h>
#include<float.h>
#include<R_ext/BLAS.h>

void computeLoop (double *Y, double *X, double *Z_delta, double *P_L, double *BF, double *sqrt_det_R,
		double *beta, double *gamma, double *phi2, double *sig2, 
		int *n, int *p, int *s, int *s_delta, int *pos_delta, int *s_L, int *pos_L, 
		int *nsim_gs, int *sum_delta) {
	char notrans='n';
	double *tmp1, *tmp2, *gamma_delta, *gamma_L;
	double alpha=1.0, beta1=0.0, f1, f2, f3, f4;
	int inc=1, i, j;

	tmp1 = (double *) calloc(*n, sizeof(double));
	tmp2 = (double *) calloc(*n, sizeof(double));
	gamma_delta = (double *) calloc(*s_delta, sizeof(double));
	gamma_L = (double *) calloc(*s_L, sizeof(double));

        /* Loops through the simulated values from the Gibbs sampler to compute the estimated 
           Bayes Factor */
	for(i=0; i<*nsim_gs; i++) {
                /* extract the gamma parameters for the right simulation iteration */
		for(j=0;j<*s_delta;j++) gamma_delta[j] = *(gamma + (*s)*i + pos_delta[j] - 1);
		for(j=0;j<*s_L;j++) gamma_L[j] = *(gamma + (*s) * i + pos_L[j] - 1);

		f1 = R_pow(phi2[i], 0.5 * (*s_L));
		if (f1 == 0.0) {
			f1 = DBL_MIN;
			printf("Warning: f1 reset to DBL_MIN\n");
 		}

		F77_CALL(dcopy)(n, Y, &inc, tmp1, &inc);
		alpha = -1.0; 
		beta1 = 1.0;

		F77_CALL(dgemv)(&notrans, n, p, &alpha, X, n, beta + (*p)*i, &inc, &beta1, tmp1, &inc);

		if(*sum_delta>1) {
			F77_CALL(dgemv)(&notrans, n, s_delta, &alpha, Z_delta, n, gamma_delta, &inc, &beta1, tmp1, &inc);
		}
		else {
			alpha = -1.0 * (*gamma_delta);
			F77_CALL(daxpy)(n, &alpha, Z_delta, &inc, tmp1, &inc);
		}

		alpha=1.0;
		beta1=0.0;
		F77_CALL(dgemv)(&notrans, n, n, &alpha, P_L, n, tmp1, &inc, &beta1, tmp2, &inc);
		f4 = F77_CALL(ddot)(n, tmp2, &inc, tmp1, &inc);
		f2 = exp((double) -0.5 * (f4/sig2[i]));
		if (f2 == 0.0) {
			f2 = DBL_MIN;
			printf("Warning: f2 reset to DBL_MIN\n");
 		}

		BF[i] = f1 * (*sqrt_det_R) * f2;
		if (BF[i] == 0.0) {
			BF[i] = DBL_MIN;
			printf("Warning: BF[i] reset to DBL_MIN\n");
 		}
	}
	

	free(tmp1);
	free(tmp2);
	free(gamma_delta);
	free(gamma_L);
}
