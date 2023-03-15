#pragma once

//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

//This project
#include "BiologicalWeightingFunction.h"
#include "NonLinearLeastSquaresUtilities.h"


class AlphaBeta_Fitter
{
	public:

		//Functions to set up the fitter
		void SetCellStudyParameters(CellStudyBWFFittingParameters survivalParams) {_survivalParams=survivalParams; _paramsSet = true; };

		//Public interface to start the fit
		double* Fit(double* initialGuess, bool weightedFit=false); 

	private:

		//Internal functions to conduct the fitting
		static int LQModel(const gsl_vector* x, void* data, gsl_vector* f);
		static double* LQModelFitting(CellStudyBWFFittingParameters& survivalParams, double* initialGuess, bool weightedFit);

		//The fitter owns its survival parameters and fitting function
		CellStudyBWFFittingParameters _survivalParams{};

		//Flags
		bool _paramsSet{false};
};
