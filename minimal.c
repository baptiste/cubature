/* compiled with the following */
/* g++ -DHCUBATURE -o minimal hcubature.c minimal.c -lm -I /usr/local/include -O2  -larmadillo -framework Accelerate */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cubature.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

/* integrand using Armadillo data structures */
colvec fa(const int ndim, const double sigma, double x)
{
  colvec results(ndim);
  int ii;
  for (ii=0; ii<ndim; ii++)
    {
      results(ii) = exp( - arma::datum::pi*sigma * ii * x*x );
    }
  return results;
}


/* wrapper of integrand for integration */
int fwrap(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    double sigma = *((double *) fdata);

    colvec res(fval, fdim, false);
    res = fa(ndim, sigma, *x);

    /* fail to access the results directly */
    // fval = res.memptr();

    // colvec res = fa(ndim, sigma, *x);
    //  int ii;
    // for (ii=0; ii<ndim; ii++) fval[ii] = res[ii];

    return 0;
}

int main(int argc, char** argv){

  const int ndim = 3;
  const int fdim = 3;

  // initialise the vectors to store integration results
  std::vector<double> integral(fdim);
  std::vector<double> error(fdim);
  double* integral_pt = &integral[0];
  double* error_pt = &error[0];
  
  double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5;

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

   hcubature(fdim, fwrap, &sigma, ndim, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, integral_pt, error_pt);

  // initialise an Armadillo matrix to use external memory
  mat  result(integral_pt, 1, fdim, false);
  // do things with the matrix...
  result.print("result:");
  
   return 0;
   
}

