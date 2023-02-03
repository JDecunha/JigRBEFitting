#include <iostream>

//This project
#include "NonLinearLeastSquaresUtilities.h"

//GNU Scientific Library
#include <gsl/gsl_blas.h>

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stdout, "iter %2zu:", iter);
  for (int i = 0; i < *(int*)params; ++i)
  {
    fprintf(stdout, " C_%i = %.4f",i,gsl_vector_get(x,i));
  }

  fprintf(stdout, ". cond(J) = %8.4f, |f(x)| = %.4f\n", 1.0 / rcond, gsl_blas_dnrm2(f));
}