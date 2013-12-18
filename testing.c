#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cubature.h"

int count = 0;
unsigned integrand_fdim = 0;
int *which_integrand = NULL;
const double radius = 0.50124145262344534123412; /* random */

/* Simple constant function */
double
fconst (double x[], size_t dim, void *params)
{
  return 1;
}

/*** f0, f1, f2, and f3 are test functions from the Monte-Carlo
     integration routines in GSL 1.6 (monte/test.c).  Copyright (c)
     1996-2000 Michael Booth, GNU GPL. ****/

/* Simple product function */
double f0 (unsigned dim, const double *x, void *params)
{
     double prod = 1.0;
     unsigned int i;
     for (i = 0; i < dim; ++i)
	  prod *= 2.0 * x[i];
     return prod;
}

#define K_2_SQRTPI 1.12837916709551257390

/* Gaussian centered at 1/2. */
double f1 (unsigned dim, const double *x, void *params)
{
     double a = *(double *)params;
     double sum = 0.;
     unsigned int i;
     for (i = 0; i < dim; i++) {
	  double dx = x[i] - 0.5;
	  sum += dx * dx;
     }
     return (pow (K_2_SQRTPI / (2. * a), (double) dim) *
	     exp (-sum / (a * a)));
}

/* double gaussian */
double f2 (unsigned dim, const double *x, void *params)
{
     double a = *(double *)params;
     double sum1 = 0.;
     double sum2 = 0.;
     unsigned int i;
     for (i = 0; i < dim; i++) {
	  double dx1 = x[i] - 1. / 3.;
	  double dx2 = x[i] - 2. / 3.;
	  sum1 += dx1 * dx1;
	  sum2 += dx2 * dx2;
     }
     return 0.5 * pow (K_2_SQRTPI / (2. * a), dim) 
	  * (exp (-sum1 / (a * a)) + exp (-sum2 / (a * a)));
}


#define K_PI 3.14159265358979323846

#include <ctype.h>
int main(int argc, char **argv)
{
     double *xmin, *xmax;
     double tol, *val, *err;
     unsigned i, dim, maxEval;

     if (argc <= 1) {
	  fprintf(stderr, "Usage: %s [dim] [reltol] [integrand] [maxeval]\n",
		  argv[0]);
	  return EXIT_FAILURE;
     }

     dim = argc > 1 ? atoi(argv[1]) : 2;
     tol = argc > 2 ? atof(argv[2]) : 1e-2;
     maxEval = argc > 4 ? atoi(argv[4]) : 0;
     
     /* parse: e.g. "x/y/z" is treated as fdim = 3, which_integrand={x,y,z} */
     if (argc <= 3) {
	  integrand_fdim = 1;
	  which_integrand = (int *) malloc(sizeof(int) * integrand_fdim);
	  which_integrand[0] = 0; /* default */
     }
     else {
	  unsigned j = 0;
	  integrand_fdim = 1;
	  for (i = 0; argv[3][i]; ++i) if (argv[3][i] == '/') ++integrand_fdim;
	  if (!integrand_fdim) {
	       fprintf(stderr, "invalid which_integrand \"%s\"", argv[3]);
	       return EXIT_FAILURE;
	  }
	  which_integrand = (int *) malloc(sizeof(int) * integrand_fdim);
	  which_integrand[0] = 0;
	  for (i = 0; argv[3][i]; ++i) {
	       if (argv[3][i] == '/')
		    which_integrand[++j] = 0;
	       else if (isdigit(argv[3][i]))
		    which_integrand[j] = 
			 which_integrand[j]*10 + argv[3][i] - '0';
	       else {
		    fprintf(stderr, "invalid which_integrand \"%s\"", argv[3]);
		    return EXIT_FAILURE;
	       }
	  }
     }
     val = (double *) malloc(sizeof(double) * integrand_fdim);
     err = (double *) malloc(sizeof(double) * integrand_fdim);

     xmin = (double *) malloc(dim * sizeof(double));
     xmax = (double *) malloc(dim * sizeof(double));
     for (i = 0; i < dim; ++i) {
	  xmin[i] = 0;
	  xmax[i] = 1;
     }

     printf("%u-dim integral, tolerance = %g\n", dim, tol);
     cubature(integrand_fdim, f_test, NULL, 
	      dim, xmin, xmax, 
	      maxEval, 0, tol, ERROR_INDIVIDUAL, val, err);
     for (i = 0; i < integrand_fdim; ++i) {
	  printf("integrand %d: integral = %0.11g, est err = %g, true err = %g\n", 
		 which_integrand[i], val[i], err[i], 
		 fabs(val[i] - exact_integral(which_integrand[i], dim, xmax)));
     }
     printf("#evals = %d\n", count);

     free(xmax);
     free(xmin);
     free(err);
     free(val);
     free(which_integrand);

     return EXIT_SUCCESS;
}
