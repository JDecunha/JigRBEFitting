#pragma once

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

class BiologicalWeightingFunction;

void H460Fitting();
void GeneralizedH460BWFFitting(BiologicalWeightingFunction fitFunc);