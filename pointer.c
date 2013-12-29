/* compiled with the following */
/* g++ -DHCUBATURE -o pointer hcubature.c pointer.c -lm -I /usr/local/include -O2  -larmadillo -framework Accelerate */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cubature.h"

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

/* vector-valued integrand using Armadillo data structures */
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

    colvec res = fa(ndim, sigma, *x);

    /* fail to access the results directly */
    /* fval = &res[0]; */  
    //typedef std::vector<double> stdvec;
    // *fval = &conv_to< stdvec >::from(res);
     int ii; 
     for (ii=0; ii<ndim; ii++) fval[ii] = res[ii]; 

    return 0;
}

int main(int argc, char** argv){

  const int fdim = 2;
  std::vector<double> memory_block(fdim);
  double* data_pointer = &memory_block[0];

  // initialise an Armadillo matrix to use external memory
  Mat<double>  result(data_pointer, 1, fdim, false);
  
  double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5;
  double *val, *err;

   err =  (double *) malloc(sizeof(double) * 2);

   /* int hcubature(unsigned fdim, integrand f, void *fdata, */
   /*               unsigned dim, const double *xmin, const double *xmax,  */
   /*               size_t maxEval, double reqAbsError, double reqRelError,  */
   /*               error_norm norm, */
   /*               double *val, double *err); */

   hcubature(fdim, fwrap, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, data_pointer, err);

   result.print("result:");
   
   return 0;
   
}

