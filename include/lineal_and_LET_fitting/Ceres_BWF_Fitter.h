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
#include "BWF_Fitting_Results.h"


class Ceres_BWF_Fitter
{
	public:

		//Functions to set up the fitter
		void SetAlphaWeightingFunction(BiologicalWeightingFunction fitFunc) {_alphaFitFunc = fitFunc; _alphaFunctionSet = true; };
		void SetBetaWeightingFunction(BiologicalWeightingFunction fitFunc) {_betaFitFunc = fitFunc; _betaFunctionSet = true; };
		void SetCellStudyParameters(CellStudyBWFFittingParameters survivalParams) {_survivalParams=survivalParams; _paramsSet = true; };

		//Functions to set constraints on fitting parameters in the alpha and beta functions
		void SetAlphaParameterLowerConstraint(int index, double constraint);
		void SetBetaParameterLowerConstraint(int index, double constraint);
		void SetAlphaParameterUpperConstraint(int index, double constraint);
		void SetPositiveConstrained(bool constrained, double penaltyWeight = 1.) {_positiveConstrained = constrained; _penaltyWeight = penaltyWeight;};


		//Function to set ceres::Solver::Options with a user defined version rather than default
		void SetFittingOptions(ceres::Solver::Options const& options) {_options=options; _optionsSet = true; };

		//These functions became necessary when we started passing constraints to the fitter
		//We need to initialize the problem before setting the constraints
		void Initialize() 
		{ 
			if(!_paramsSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting survival params. Failure." << std::endl; throw; }
			if(!_alphaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  alpha biological weighting function to fit. Failure." << std::endl; throw; }
			if(!_betaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  beta biological weighting function to fit. Failure." << std::endl; throw; }
			_problem = ceres::Problem(); SetupResidualBlocks(); _initialized = true; _alreadyRun = false; 
		};

		//Public interface to start the fit
		BWF_Fitting_Results Fit(); 

		//Static helper functions
		static void CheckFunctionNegativity(BiologicalWeightingFunction const& BWF, double lower = 0.1, double upper = 100); 

	private:

		struct Generalized_BWF_Residual
		{
			Generalized_BWF_Residual(double dose, double SF, TH1D const& linealSpectrum, BiologicalWeightingFunction const& alphaFunc, BiologicalWeightingFunction const& betaFunc) 
				: _dose(dose), _SF(SF), _linealSpectrum(linealSpectrum), _alphaFunc(alphaFunc), _betaFunc(betaFunc) 
			{

			}

			bool operator()(double const* const* parameters, double* residual) const 
			{

				double alphaPredicted = 0.; double betaPredicted = 0.;

				for(int i = 1; i <= _linealSpectrum.GetNbinsX(); ++i) //Iterate over the lineal energy spectrum
				{
					//Get the spectrum value
					double width = _linealSpectrum.GetBinWidth(i);
					double center = _linealSpectrum.GetBinCenter(i);
					double value = _linealSpectrum.GetBinContent(i);

					//Accumulate the predicted value of alpha as we integrate the function
					double ryAlphaVal = _alphaFunc.GetValue(parameters[0],center);
					alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

					double ryBetaVal = _betaFunc.GetValue(parameters[1], center);
					betaPredicted += ryBetaVal*value*width;
				}

				//Logarithmic residual
				double survivalPredicted = -( (alphaPredicted*_dose) + (betaPredicted*_dose*_dose) );
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
				TH1D const& _linealSpectrum;
				BiologicalWeightingFunction const _alphaFunc;
				BiologicalWeightingFunction const _betaFunc;
		};

		struct BWF_Negative_Penalty
		{
			BWF_Negative_Penalty(BiologicalWeightingFunction const& alphaFunc, BiologicalWeightingFunction const& betaFunc, double penaltyWeight = 1) 
				: _alphaFunc(alphaFunc), _betaFunc(betaFunc), _penaltyWeight(penaltyWeight)
			{

			}

			bool operator()(double const* const* parameters, double* residual) const 
			{
				double penaltyTerm = 0.;

				for(double center = 0.1; center <= 120; center+=0.1) //Iterate over the lineal energy spectrum
				{
					//Accumulate the predicted value of alpha as we integrate the function
					double ryAlphaVal = _alphaFunc.GetValue(parameters[0],center);
					double ryBetaVal = _betaFunc.GetValue(parameters[1], center);

					//Evaluate penalty terms
					if(ryAlphaVal < 0)
					{
						penaltyTerm += -ryAlphaVal*_penaltyWeight;
					}
					//Evaluate penalty terms
					if(ryBetaVal < 0)
					{
						penaltyTerm += -ryBetaVal*_penaltyWeight;
					}
				}

				//The residuals are automatically squared in the cost function. So I'm taking SQRT to make it linear.
				residual[0] = std::sqrt(penaltyTerm)*std::sqrt(2);
		    	return true;
		 	}

			private:
				//Everything is const so once a residual block is defined it can't be changed
				double const _penaltyWeight;
				BiologicalWeightingFunction const _alphaFunc;
				BiologicalWeightingFunction const _betaFunc;
		};


		//Internal functions to conduct the fitting
		void GeneralizedBWFFitting();
		void SetupResidualBlocks();

		//The fitter owns its survival parameters and fitting function
		CellStudyBWFFittingParameters _survivalParams{};
		BiologicalWeightingFunction _alphaFitFunc{};
		BiologicalWeightingFunction _betaFitFunc{};

		//The fitter owns the problem it solves, the options associated with it, and the results of the last fit
		ceres::Problem _problem{};
		ceres::Solver::Options _options;
		BWF_Fitting_Results results{};

		//Flags
		bool _paramsSet{false};
		bool _alphaFunctionSet{false};
		bool _betaFunctionSet{false};
		bool _positiveConstrained{false};
		bool _optionsSet{false};
		bool _initialized{false};
		bool _alreadyRun{false};
		double _penaltyWeight{1.};
};



class Ceres_BWF_Fitting_Callback : public ceres::IterationCallback 
{
	public:
		explicit Ceres_BWF_Fitting_Callback(BiologicalWeightingFunction* alphaFunc, BiologicalWeightingFunction* betaFunc) 
		: _pAlphaFunc(alphaFunc), _pBetaFunc(betaFunc)
		{}
		
		~Ceres_BWF_Fitting_Callback() {}

		ceres::CallbackReturnType operator() (const ceres::IterationSummary& summary) 
		{
			std::cout << "\u001b[31m\e[1mIteration #: \e[0m\u001b[0m" << summary.iteration;
			std::cout << "\u001b[34m\e[1m Alpha: \e[0m\u001b[0m"; _pAlphaFunc->PrintFitParams();
			std::cout << "\u001b[34m\e[1m  Beta: \e[0m\u001b[0m"; _pBetaFunc->PrintFitParams();
			std::cout << std::endl;

			return ceres::SOLVER_CONTINUE;
		}

	private:
		BiologicalWeightingFunction* _pAlphaFunc;
		BiologicalWeightingFunction* _pBetaFunc;
};
