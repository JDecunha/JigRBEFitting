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


class Ceres_BWF_Fitter
{
	public:

		//Functions to set up the fitter
		void SetAlphaWeightingFunction(BiologicalWeightingFunction fitFunc) {_alphaFitFunc = fitFunc; _alphaFunctionSet = true; };
		void SetBetaWeightingFunction(BiologicalWeightingFunction fitFunc) {_betaFitFunc = fitFunc; _betaFunctionSet = true; };
		void SetCellStudyParameters(CellStudyBWFFittingParameters survivalParams) {_survivalParams=survivalParams; _paramsSet = true; };

		//Public interface to start the fit
		void Fit(); 

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

				//calculate SF
				double survivalPredicted = ( (alphaPredicted*_dose) + (betaPredicted*_dose*_dose) );
				survivalPredicted = std::exp(-survivalPredicted);

				//Return the residual
				residual[0] = _SF - survivalPredicted;
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

		//Internal functions to conduct the fitting
		void GeneralizedBWFFitting();
		void SetupResidualBlocks();

		//The fitter owns its survival parameters and fitting function
		CellStudyBWFFittingParameters _survivalParams{};
		BiologicalWeightingFunction _alphaFitFunc{};
		BiologicalWeightingFunction _betaFitFunc{};
		ceres::Problem _problem{};

		//Flags
		bool _paramsSet{false};
		bool _alphaFunctionSet{false};
		bool _betaFunctionSet{false};
};
