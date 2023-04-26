#pragma once

//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

//Ceres fitting tools
#include "ceres/ceres.h"

//This project
#include "BiologicalWeightingFunction.h"
#include "NonLinearLeastSquaresUtilities.h"
#include "Survival_Fitting_Results.h"


class Ceres_Survival_Fitter
{
	public:

		//Functions to set up the fitter

		// void SetAlphaWeightingFunction(BiologicalWeightingFunction fitFunc) {_alphaFitFunc = fitFunc; _alphaFunctionSet = true; };
		// void SetBetaWeightingFunction(BiologicalWeightingFunction fitFunc) {_betaFitFunc = fitFunc; _betaFunctionSet = true; };
		void SetCellStudyParameters(CellStudyBWFFittingParameters survivalParams) {_survivalParams=survivalParams; _paramsSet = true; };

		//Functions to set constraints on fitting parameters in the alpha and beta functions
		// void SetAlphaParameterLowerConstraint(int index, double constraint);
		// void SetBetaParameterLowerConstraint(int index, double constraint);
		// void SetAlphaParameterUpperConstraint(int index, double constraint);
		// void SetPositiveConstrained(bool constrained, double penaltyWeight = 1.) {_positiveConstrained = constrained; _penaltyWeight = penaltyWeight;};


		//Function to set ceres::Solver::Options with a user defined version rather than default
		void SetFittingOptions(ceres::Solver::Options const& options) {_options=options; _optionsSet = true; };

		//These functions became necessary when we started passing constraints to the fitter
		//We need to initialize the problem before setting the constraints
		void Initialize() 
		{ 
			if(!_paramsSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting survival params. Failure." << std::endl; throw; }
			// if(!_alphaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  alpha biological weighting function to fit. Failure." << std::endl; throw; }
			// if(!_betaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  beta biological weighting function to fit. Failure." << std::endl; throw; }
			_problem = ceres::Problem(); SetupResidualBlocks(); _initialized = true; _alreadyRun = false; 
		};

		//Public interface to start the fit
		Survival_Fitting_Results Fit(); 

		//Static helper functions
		// static void CheckFunctionNegativity(BiologicalWeightingFunction const& BWF, double lower = 0.1, double upper = 100); 

	private:

		struct Generalized_Survival_Residual
		{
			Generalized_Survival_Residual(double dose, double SF) 
				: _dose(dose), _SF(SF)
			{}

			bool operator()(double const* const* parameters, double* residual) const 
			{

				//Logarithmic residual
				double alpha = parameters[0][0]; double beta = parameters[1][0];
				double survivalPredicted = -( (alpha*_dose) + (beta*_dose*_dose) );
				residual[0] = std::log(_SF) - survivalPredicted;

				//Simple surviving fraction residual
				// double survivalPredicted = -( (alphaPredicted*_dose) + (betaPredicted*_dose*_dose) );
				// survivalPredicted = std::exp(survivalPredicted);
				// residual[0] = _SF - survivalPredicted;

				//Percent difference residual
				// double survivalPredicted = -( (alphaPredicted*_dose) + (betaPredicted*_dose*_dose) );
				// survivalPredicted = std::exp(survivalPredicted);
				// residual[0] = (survivalPredicted - _SF)/survivalPredicted;

		    	return true;
		 	}

			private:

				//Everything is const so once a residual block is defined it can't be changed
				double const _dose;
				double const _SF;
		};


		//Internal functions to conduct the fitting
		void GeneralizedSurvivalFitting();
		void SetupResidualBlocks();

		//The fitter owns its survival parameters and fitting function
		CellStudyBWFFittingParameters _survivalParams{};

		//The fitter owns the problem it solves, the options associated with it, and the results of the last fit
		ceres::Problem _problem{};
		ceres::Solver::Options _options;
		Survival_Fitting_Results results;

		//Flags
		bool _paramsSet{false};
		bool _optionsSet{false};
		bool _initialized{false};
		bool _alreadyRun{false};
		double _penaltyWeight{1.};
};



class Ceres_Survival_Fitting_Callback : public ceres::IterationCallback 
{
	public:
		explicit Ceres_Survival_Fitting_Callback() 	{}
		
		~Ceres_Survival_Fitting_Callback() {}

		ceres::CallbackReturnType operator() (const ceres::IterationSummary& summary) 
		{
			std::cout << "\u001b[31m\e[1mIteration #: \e[0m\u001b[0m" << summary.iteration;
			std::cout << std::endl;

			return ceres::SOLVER_CONTINUE;
		}
};
