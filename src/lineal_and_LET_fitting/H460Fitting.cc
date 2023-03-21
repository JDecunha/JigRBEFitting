//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

//ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"

//This project
#include "ProtonSpectra.h" //Proton KESpectra
#include "LinealSpectra.h" //Lineal energy spectra from KE Spectra
//Functions for fitting and plotting
#include "SurvivalPlotting.h"
// #include "BWF_Fitter.h"
// #include "BWF_Fitter_Beta.h"
#include "LET_Fitter.h"
#include "BWF_Fitter_AlphaBeta.h"
#include "AlphaBeta_Fitter.h"
#include "BWF_Plotter.h"
//This file
#include "H460Fitting.h"

//I put the survival parameters at global scope
// CellStudySurvivalParameters H460Params{};
//The version below works with having a BWF for beta too.
extern CellStudyBWFFittingParameters H460FittingParams;

// void BetaFromPhotonFitting()
// {
// 	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
// 	SetupH460SurvivalParameters();

// 	//Plot
// 	//Setup the canvas
// 	gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	TCanvas* c = new TCanvas("c","c");
// 	c->SetCanvasSize(9000, 5000);
// 	c->SetFillStyle(4000);
// 	c->SetFrameFillStyle(4000);
// 	c->Divide(4,3,0.000000005,0.001);
// 	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
// 	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
// 	legend->SetTextSize(0.05);

// 	//Setup the marker attributes
// 	TAttLine lineStyle{};
// 	lineStyle.SetLineColor(kGreen+2);
// 	lineStyle.SetLineWidth(3);
// 	lineStyle.SetLineStyle(2);

// 	//Plot the survival data
// 	SurvivalDataMultigraph(c, legend, H460Params);

// 	//
// 	//Set up the BWFs
// 	//

// 	//Set up the fitter
// 	BWF_Fitter_AlphaBeta fitter{};
// 	fitter.SetCellStudyParameters(H460FittingParamsNewAndImproved);

// 	//Create photon beta "BWF"
// 	BiologicalWeightingFunction photonBeta;
// 	photonBeta.SetWeightingFunction([](double* params, double linealEnergy) {return 0.097;}, 0);

// 	//Create fixed "BWF"
// 	BiologicalWeightingFunction FixedBWF;
// 	FixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[0];}, 1);

// 	//Create linear BWF
// 	BiologicalWeightingFunction LinearandFixedBWF;
// 	LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

// 	//Exponential growth BWF
// 	BiologicalWeightingFunction ExponentialGrowthBWF;
// 	ExponentialGrowthBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

// 	//Decaying to zero BWF
// 	BiologicalWeightingFunction ExponentialToZeroBWF;
// 	ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

// 	BiologicalWeightingFunction ExponentialToZeroRestrictedBWF;
// 	ExponentialToZeroRestrictedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return 0.15*(1-std::exp(-params[0]/linealEnergy)) ;}, 1);

// 	// //Decaying to zero BWF, simple
// 	// BiologicalWeightingFunction ExponentialToZeroBWF;
// 	// ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(params[0]/linealEnergy)) ;}, 2);

// 	//Simple exponential decay BWF
// 	BiologicalWeightingFunction ExponentialDecayBWF;
// 	ExponentialDecayBWF.SetWeightingFunction([](double* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

// 	//Create quadratic BWF
// 	BiologicalWeightingFunction QuadraticLinearFixedBWF;
// 	QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearFixedBWF;
// 	CubicQuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearBWF;
// 	CubicQuadraticLinearBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy*linealEnergy)+(params[1]*linealEnergy*linealEnergy)+(params[0]*linealEnergy));}, 3);

// 	//Create factored quadratic BWF
// 	BiologicalWeightingFunction QuadraticFactoredBWF;
// 	QuadraticFactoredBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[1]*(linealEnergy-params[2])*(linealEnergy-params[2]))+params[0]);}, 3);

// 	//BiologicalWeightingFunction GaussianBWF;
// 	BiologicalWeightingFunction GaussianBWF;
// 	GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
// 	{ 
// 		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
// 	}
// 	, 2);

// 	///
// 	/// The actual fitting
// 	///

// 	// Linear fitting
// 	fitter.SetAlphaWeightingFunction(LinearandFixedBWF);
// 	fitter.SetBetaWeightingFunction(photonBeta);
// 	double linearandBetaFixedInitialGuess [] = {0.9, 0.2};
// 	double* LinearandBetaFixedParams = fitter.Fit(linearandBetaFixedInitialGuess,true);
// 	//Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kBlue+2);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Linear", H460FittingParamsNewAndImproved, LinearandFixedBWF, photonBeta, LinearandBetaFixedParams);

// 	//Quadratic fitting
// 	fitter.SetAlphaWeightingFunction(QuadraticLinearFixedBWF);
// 	fitter.SetBetaWeightingFunction(photonBeta);
// 	double quadraticAndBetaFixedInitialGuess [] = {0.2, 0., 0.01};
// 	double* quadraticAndBetaFixedParams = fitter.Fit(quadraticAndBetaFixedInitialGuess,true);
// 	// Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kRed+2);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, photonBeta, quadraticAndBetaFixedParams);

// 	//Cubic Fitting
// 	fitter.SetAlphaWeightingFunction(CubicQuadraticLinearFixedBWF);
// 	fitter.SetBetaWeightingFunction(photonBeta);	
// 	double cubicQuadraticAndBetaFixedInitialGuess [] = {0.06, 0.1, -0.063, 0.002};
// 	double* cubicQuadraticAndBetaFixedParams = fitter.Fit(cubicQuadraticAndBetaFixedInitialGuess,true);
// 	// // Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kGreen+4);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic", H460FittingParamsNewAndImproved, CubicQuadraticLinearFixedBWF, QuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams);

// 	// //Cubic, no linear offset Fitting
// 	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearBWF);
// 	// fitter.SetBetaWeightingFunction(FixedBWF);	
// 	// double cubicQuadraticLinearAndBetaFixedInitialGuess [] = {0.1177, -0.02, 0.003, 0.13};
// 	// double* cubicQuadraticLinearAndBetaFixedParams = fitter.Fit(cubicQuadraticLinearAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kBlue+4);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic, sans offset", H460FittingParamsNewAndImproved, CubicQuadraticLinearBWF, FixedBWF, cubicQuadraticLinearAndBetaFixedParams);

// 	// //Gaussian Fitting
// 	fitter.SetAlphaWeightingFunction(GaussianBWF);
// 	fitter.SetBetaWeightingFunction(photonBeta);	
// 	double gaussianInitialGuess [] = {0.0008, 95.4568}; //5, 20, 0.13
// 	double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);


// 	// //Save
// 	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_Fits_Linear_and_quadratic_H460.jpg";
// 	// c->SaveAs((TString)outputName); 


// 	// //Plotting the BWFs themselves
// 	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	// TCanvas* c2 = new TCanvas("c","c");
// 	// c2->SetCanvasSize(4000, 2400);
// 	// // auto legend2 = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
// 	// auto legend2 = new TLegend(0.72,0.14,0.72+0.2,0.14+0.16);//x start, y start, x end, yend
// 	// legend2->SetTextSize(0.04);

// 	// //Setup the marker attributes
// 	// TAttLine lineStyle2{};
	
// 	// lineStyle2.SetLineWidth(5);
// 	// lineStyle2.SetLineStyle(1);

// 	// // //Plot the BWFs
// 	// // lineStyle2.SetLineColor(kGreen+4);
// 	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic", CubicQuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams, "AL", 0., 100.);
// 	// // // lineStyle2.SetLineColor(kBlue+4);
// 	// // // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic, sans offset", CubicQuadraticLinearBWF, cubicQuadraticLinearAndBetaFixedParams, "L", 0., 100.);
// 	// lineStyle2.SetLineColor(kBlue+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Linear", LinearandFixedBWF, LinearandBetaFixedParams,"AL");
// 	// lineStyle2.SetLineColor(kRed+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Quadratic", QuadraticLinearFixedBWF, quadraticAndBetaFixedParams, "L");

// 	// //Save
// 	// std::string outputName2 = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_BWFs_Linear_and_Quadratic_H460.jpg";
// 	// c2->SaveAs((TString)outputName2); 
// }

// void LETFixedBetaFitting()
// {
// 	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
// 	SetupH460SurvivalParameters();

// 	//Plot
// 	//Setup the canvas
// 	gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	TCanvas* c = new TCanvas("c","c");
// 	c->SetCanvasSize(9000, 5000);
// 	c->SetFillStyle(4000);
// 	c->SetFrameFillStyle(4000);
// 	c->Divide(4,3,0.000000005,0.001);
// 	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
// 	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
// 	legend->SetTextSize(0.05);

// 	// //Setup the marker attributes
// 	TAttLine lineStyle{};
// 	lineStyle.SetLineColor(kGreen+2);
// 	lineStyle.SetLineWidth(3);
// 	lineStyle.SetLineStyle(2);

// 	//Plot the survival data
// 	SurvivalDataMultigraph(c, legend, H460Params);

// 	//
// 	//Set up the BWFs
// 	//

// 	//Set up the fitter
// 	LET_Fitter fitter{};
// 	fitter.SetCellStudyParameters(H460FittingParamsNewAndImproved);

// 	//Create fixed "BWF"
// 	BiologicalWeightingFunction FixedBWF;
// 	FixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[0];}, 1);

// 	//Create linear BWF
// 	BiologicalWeightingFunction LinearandFixedBWF;
// 	LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

// 	//Exponential growth BWF
// 	BiologicalWeightingFunction ExponentialGrowthBWF;
// 	ExponentialGrowthBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

// 	//Decaying to zero BWF
// 	BiologicalWeightingFunction ExponentialToZeroBWF;
// 	ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

// 	BiologicalWeightingFunction ExponentialToZeroRestrictedBWF;
// 	ExponentialToZeroRestrictedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return 0.15*(1-std::exp(-params[0]/linealEnergy)) ;}, 1);

// 	// //Decaying to zero BWF, simple
// 	// BiologicalWeightingFunction ExponentialToZeroBWF;
// 	// ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(params[0]/linealEnergy)) ;}, 2);

// 	//Simple exponential decay BWF
// 	BiologicalWeightingFunction ExponentialDecayBWF;
// 	ExponentialDecayBWF.SetWeightingFunction([](double* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

// 	//Create quadratic BWF
// 	BiologicalWeightingFunction QuadraticLinearFixedBWF;
// 	QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearFixedBWF;
// 	CubicQuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearBWF;
// 	CubicQuadraticLinearBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy*linealEnergy)+(params[1]*linealEnergy*linealEnergy)+(params[0]*linealEnergy));}, 3);

// 	//Create factored quadratic BWF
// 	BiologicalWeightingFunction QuadraticFactoredBWF;
// 	QuadraticFactoredBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[1]*(linealEnergy-params[2])*(linealEnergy-params[2]))+params[0]);}, 3);

// 	//BiologicalWeightingFunction GaussianBWF;
// 	BiologicalWeightingFunction GaussianBWF;
// 	GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
// 	{ 
// 		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
// 	}
// 	, 2);

// 	///
// 	/// The actual fitting
// 	///

// 	//Linear fitting
// 	fitter.SetAlphaFunction(LinearandFixedBWF);
// 	fitter.SetBetaFunction(FixedBWF);
// 	double linearandBetaFixedInitialGuess [] = {0.9, 0.2, 0.09};
// 	double* LinearandBetaFixedParams = fitter.Fit(linearandBetaFixedInitialGuess,true);
// 	//Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kBlue+2);
// 	GeneralizedLETMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Linear", H460FittingParamsNewAndImproved, LinearandFixedBWF, FixedBWF, LinearandBetaFixedParams);

// 	// //Quadratic fitting
// 	fitter.SetAlphaFunction(QuadraticLinearFixedBWF);
// 	fitter.SetBetaFunction(FixedBWF);
// 	double quadraticAndBetaFixedInitialGuess [] = {0.2, 0., 0.01, 0.12};
// 	double* quadraticAndBetaFixedParams = fitter.Fit(quadraticAndBetaFixedInitialGuess,true);
// 	// Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kRed+2);
// 	GeneralizedLETMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, FixedBWF, quadraticAndBetaFixedParams);

// 	//Cubic Fitting
// 	fitter.SetAlphaFunction(CubicQuadraticLinearFixedBWF);
// 	fitter.SetBetaFunction(FixedBWF);	
// 	double cubicQuadraticAndBetaFixedInitialGuess [] = {0.06, 0.1, -0.063, 0.002, 0.13};
// 	double* cubicQuadraticAndBetaFixedParams = fitter.Fit(cubicQuadraticAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kGreen+4);
// 	GeneralizedLETMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic", H460FittingParamsNewAndImproved, CubicQuadraticLinearFixedBWF, FixedBWF, cubicQuadraticAndBetaFixedParams);

// 	// //Cubic, no linear offset Fitting
// 	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearBWF);
// 	// fitter.SetBetaWeightingFunction(FixedBWF);	
// 	// double cubicQuadraticLinearAndBetaFixedInitialGuess [] = {0.1177, -0.02, 0.003, 0.13};
// 	// double* cubicQuadraticLinearAndBetaFixedParams = fitter.Fit(cubicQuadraticLinearAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kBlue+4);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic, sans offset", H460FittingParamsNewAndImproved, CubicQuadraticLinearBWF, FixedBWF, cubicQuadraticLinearAndBetaFixedParams);

// 	// //Gaussian Fitting
// 	// fitter.SetAlphaWeightingFunction(GaussianBWF);
// 	// fitter.SetBetaWeightingFunction(FixedBWF);	
// 	// double gaussianInitialGuess [] = {0.0008, 95.4568, 0.2061}; //5, 20, 0.13
// 	// double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);


// 	// //Save
// 	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_AlphaFits_LET_polynomials_H460.jpg";
// 	c->SaveAs((TString)outputName); 


// 	// //Plotting the BWFs themselves
// 	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	// TCanvas* c2 = new TCanvas("c","c");
// 	// c2->SetCanvasSize(4000, 2400);
// 	// auto legend2 = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
// 	// // auto legend2 = new TLegend(0.72,0.14,0.72+0.2,0.14+0.16);//x start, y start, x end, yend
// 	// legend2->SetTextSize(0.04);

// 	// //Setup the marker attributes
// 	// TAttLine lineStyle2{};
	
// 	// lineStyle2.SetLineWidth(5);
// 	// lineStyle2.SetLineStyle(1);

// 	// // //Plot the BWFs
// 	// // lineStyle2.SetLineColor(kGreen+4);
// 	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic", CubicQuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams, "AL", 0., 120.);
// 	// // // lineStyle2.SetLineColor(kBlue+4);
// 	// // // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic, sans offset", CubicQuadraticLinearBWF, cubicQuadraticLinearAndBetaFixedParams, "L", 0., 100.);
// 	// lineStyle2.SetLineColor(kBlue+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Linear", LinearandFixedBWF, LinearandBetaFixedParams,"AL", 0., 120.);
// 	// lineStyle2.SetLineColor(kRed+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Quadratic", QuadraticLinearFixedBWF, quadraticAndBetaFixedParams, "L", 0., 120.);

// 	// //Save
// 	// std::string outputName2 = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/LET_AlphaFunctions_FixedBeta_nocubic_H460.jpg";
// 	// c2->SaveAs((TString)outputName2); 
// }

void FixedBetaFitting()
{
	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
	// SetupH460SurvivalParameters();

	//Plot
	//Setup the canvas
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	legend->SetTextSize(0.05);

	//Setup the marker attributes
	TAttLine lineStyle{};
	lineStyle.SetLineColor(kGreen+2);
	lineStyle.SetLineWidth(3);
	lineStyle.SetLineStyle(2);

	//Plot the survival data
	// SurvivalDataMultigraph(c, legend, H460Params);

	//
	//Set up the BWFs
	//

	//Set up the fitter
	BWF_Fitter_AlphaBeta fitter{};
	fitter.SetCellStudyParameters(H460FittingParams);

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	// //Create photon beta "BWF"
	// BiologicalWeightingFunction photonBeta;
	// photonBeta.SetWeightingFunction([](double* params, double linealEnergy) {return 0.097;}, 0);

	//Create linear BWF
	BiologicalWeightingFunction LinearandFixedBWF;
	LinearandFixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	// //Exponential growth BWF
	// BiologicalWeightingFunction ExponentialGrowthBWF;
	// ExponentialGrowthBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	// //Decaying to zero BWF
	// BiologicalWeightingFunction ExponentialToZeroBWF;
	// ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	// BiologicalWeightingFunction ExponentialToZeroRestrictedBWF;
	// ExponentialToZeroRestrictedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return 0.15*(1-std::exp(-params[0]/linealEnergy)) ;}, 1);

	// // //Decaying to zero BWF, simple
	// // BiologicalWeightingFunction ExponentialToZeroBWF;
	// // ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(params[0]/linealEnergy)) ;}, 2);

	// //Simple exponential decay BWF
	// BiologicalWeightingFunction ExponentialDecayBWF;
	// ExponentialDecayBWF.SetWeightingFunction([](double* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	// //Create quadratic BWF
	// BiologicalWeightingFunction QuadraticLinearFixedBWF;
	// QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	// //Create cubic BWF
	// BiologicalWeightingFunction CubicQuadraticLinearFixedBWF;
	// CubicQuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	// //Create cubic BWF
	// BiologicalWeightingFunction CubicQuadraticLinearBWF;
	// CubicQuadraticLinearBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy*linealEnergy)+(params[1]*linealEnergy*linealEnergy)+(params[0]*linealEnergy));}, 3);

	// //Create factored quadratic BWF
	// BiologicalWeightingFunction QuadraticFactoredBWF;
	// QuadraticFactoredBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[1]*(linealEnergy-params[2])*(linealEnergy-params[2]))+params[0]);}, 3);

	// //BiologicalWeightingFunction GaussianBWF;
	// BiologicalWeightingFunction GaussianBWF;
	// GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
	// { 
	// 	return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	// }
	// , 2);

	// //Three parameter sigmoid;
	// BiologicalWeightingFunction ThreeParameterSigmoid;
	// ThreeParameterSigmoid.SetWeightingFunction( [] (double* params, double linealEnergy) 
	// { 
	// 	return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2]))));
	// }
	// , 3);

	// BiologicalWeightingFunction MorstinBWF;
	// MorstinBWF.SetWeightingFunction([](double* params, double linealEnergy) 
	// {
	// 	return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	// }, 4);

	///
	/// The actual fitting
	///

	// //Linear fitting
	fitter.SetAlphaWeightingFunction(LinearandFixedBWF);
	fitter.SetBetaWeightingFunction(FixedBWF);
	double linearandBetaFixedInitialGuess [] = {0.9, 0.2, 0.09};
	double* LinearandBetaFixedParams = fitter.Fit(linearandBetaFixedInitialGuess,false);
	// //Plot the BWF prediction on the survival data
	lineStyle.SetLineColor(kBlue+2);
	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParamsNewAndImproved, LinearandFixedBWF, FixedBWF, LinearandBetaFixedParams);

	// //Quadratic fitting
	// fitter.SetAlphaWeightingFunction(QuadraticLinearFixedBWF);
	// fitter.SetBetaWeightingFunction(FixedBWF);
	// double quadraticAndBetaFixedInitialGuess [] = {0.2, 0., 0.01, 0.12};
	// double* quadraticAndBetaFixedParams = fitter.Fit(quadraticAndBetaFixedInitialGuess,true);
	// // Plot the BWF prediction on the survival data
	// lineStyle.SetLineColor(kRed+2);
	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, FixedBWF, quadraticAndBetaFixedParams);

	// // //Cubic Fitting
	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearFixedBWF);
	// fitter.SetBetaWeightingFunction(FixedBWF);	
	// double cubicQuadraticAndBetaFixedInitialGuess [] = {0.06, 0.1, -0.063, 0.002, 0.13};
	// double* cubicQuadraticAndBetaFixedParams = fitter.Fit(cubicQuadraticAndBetaFixedInitialGuess,true);
	// // Plot the BWF prediction on the survival data
	// lineStyle.SetLineColor(kGreen+4);
	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParamsNewAndImproved, CubicQuadraticLinearFixedBWF, FixedBWF, cubicQuadraticAndBetaFixedParams);

	// //Cubic, no linear offset Fitting
	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearBWF);
	// fitter.SetBetaWeightingFunction(FixedBWF);	
	// double cubicQuadraticLinearAndBetaFixedInitialGuess [] = {0.1177, -0.02, 0.003, 0.13};
	// double* cubicQuadraticLinearAndBetaFixedParams = fitter.Fit(cubicQuadraticLinearAndBetaFixedInitialGuess,true);
	// // // Plot the BWF prediction on the survival data
	// lineStyle.SetLineColor(kBlue+4);
	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic, sans offset", H460FittingParamsNewAndImproved, CubicQuadraticLinearBWF, FixedBWF, cubicQuadraticLinearAndBetaFixedParams);

	// //Gaussian Fitting
	// fitter.SetAlphaWeightingFunction(GaussianBWF);
	// fitter.SetBetaWeightingFunction(FixedBWF);	
	// double gaussianInitialGuess [] = {0.0008, 95.4568, 0.2061}; //5, 20, 0.13
	// double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);

	// //Sigmoid fitting
	// fitter.SetAlphaWeightingFunction(ThreeParameterSigmoid);
	// fitter.SetBetaWeightingFunction(FixedBWF);
	// double ThreeParamSigmoidGuess [] = {6740., 0.02, 400., 0.12};
	// double* ThreeParamSigmoidFitParams = fitter.Fit(ThreeParamSigmoidGuess,true);
	// // Plot the BWF prediction on the survival data
	// // lineStyle.SetLineColor(kRed+2);
	// // GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, FixedBWF, quadraticAndBetaFixedParams);

	// //Morstin fitting
	// fitter.SetAlphaWeightingFunction(MorstinBWF);
	// fitter.SetBetaWeightingFunction(photonBeta);
	// // double MorstinGuess [] = {11460.000000, 2.5*std::pow(10,-6), 2.1*std::pow(10,-5), 2.*std::pow(10,-7)};
	// double MorstinGuess [] = {11460.000000, 0., 2.1*std::pow(10,-5), 1.*std::pow(10,-7)};
	// double* MorstinFitParams = fitter.Fit(MorstinGuess,true);


	// //Save
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_Fits_H460_CubicResiduals.jpg";
	// c->SaveAs((TString)outputName); 


	//Plotting the BWFs themselves
	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(4000, 2400);
	// auto legend2 = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
	// // auto legend2 = new TLegend(0.72,0.14,0.72+0.2,0.14+0.16);//x start, y start, x end, yend
	// legend2->SetTextSize(0.04);

	// //Setup the marker attributes
	// TAttLine lineStyle2{};
	
	// lineStyle2.SetLineWidth(5);
	// lineStyle2.SetLineStyle(1);

	// // //Plot the BWFs
	// lineStyle2.SetLineColor(kGreen+4);
	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic", CubicQuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams, "AL", 0., 120.);
	// // // lineStyle2.SetLineColor(kBlue+4);
	// // // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic, sans offset", CubicQuadraticLinearBWF, cubicQuadraticLinearAndBetaFixedParams, "L", 0., 100.);
	// lineStyle2.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Linear", LinearandFixedBWF, LinearandBetaFixedParams,"L", 0., 120.);
	// lineStyle2.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Quadratic", QuadraticLinearFixedBWF, quadraticAndBetaFixedParams, "L", 0., 120.);

	// // lineStyle2.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Sigmoid-Test", ThreeParameterSigmoid, ThreeParamSigmoidGuess, "L");
	// // lineStyle2.SetLineColor(kViolet+2);
	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Sigmoid-Real", ThreeParameterSigmoid, ThreeParamSigmoidFitParams, "L");

	// // lineStyle2.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Morstin-Initial", MorstinBWF, MorstinGuess, "L");
	// // lineStyle2.SetLineColor(kViolet+2);
	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Morstin-Real", MorstinBWF, MorstinFitParams, "L");

	// //Save
	// std::string outputName2 = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_BWFs_Polynomial_H460.jpg";
	// c2->SaveAs((TString)outputName2); 
}

// void LinearBetaFitting()
// {
// 	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
// 	SetupH460SurvivalParameters();

// 	//Plot
// 	//Setup the canvas
// 	gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	TCanvas* c = new TCanvas("c","c");
// 	c->SetCanvasSize(9000, 5000);
// 	c->SetFillStyle(4000);
// 	c->SetFrameFillStyle(4000);
// 	c->Divide(4,3,0.000000005,0.001);
// 	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
// 	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
// 	legend->SetTextSize(0.05);

// 	//Setup the marker attributes
// 	TAttLine lineStyle{};
// 	lineStyle.SetLineColor(kGreen+2);
// 	lineStyle.SetLineWidth(3);
// 	lineStyle.SetLineStyle(2);

// 	//Plot the survival data
// 	SurvivalDataMultigraph(c, legend, H460Params);

// 	//
// 	//Set up the BWFs
// 	//

// 	//Set up the fitter
// 	BWF_Fitter_AlphaBeta fitter{};
// 	fitter.SetCellStudyParameters(H460FittingParamsNewAndImproved);

// 	//Create fixed "BWF"
// 	BiologicalWeightingFunction FixedBWF;
// 	FixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[0];}, 1);

// 	//Create linear BWF
// 	BiologicalWeightingFunction LinearandFixedBWF;
// 	LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

// 	//Exponential growth BWF
// 	BiologicalWeightingFunction ExponentialGrowthBWF;
// 	ExponentialGrowthBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

// 	//Decaying to zero BWF
// 	BiologicalWeightingFunction ExponentialToZeroBWF;
// 	ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

// 	BiologicalWeightingFunction ExponentialToZeroRestrictedBWF;
// 	ExponentialToZeroRestrictedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return 0.15*(1-std::exp(-params[0]/linealEnergy)) ;}, 1);

// 	// //Decaying to zero BWF, simple
// 	// BiologicalWeightingFunction ExponentialToZeroBWF;
// 	// ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(params[0]/linealEnergy)) ;}, 2);

// 	//Simple exponential decay BWF
// 	BiologicalWeightingFunction ExponentialDecayBWF;
// 	ExponentialDecayBWF.SetWeightingFunction([](double* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

// 	//Create quadratic BWF
// 	BiologicalWeightingFunction QuadraticLinearFixedBWF;
// 	QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearFixedBWF;
// 	CubicQuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearBWF;
// 	CubicQuadraticLinearBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy*linealEnergy)+(params[1]*linealEnergy*linealEnergy)+(params[0]*linealEnergy));}, 3);

// 	//Create factored quadratic BWF
// 	BiologicalWeightingFunction QuadraticFactoredBWF;
// 	QuadraticFactoredBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[1]*(linealEnergy-params[2])*(linealEnergy-params[2]))+params[0]);}, 3);

// 	//BiologicalWeightingFunction GaussianBWF;
// 	BiologicalWeightingFunction GaussianBWF;
// 	GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
// 	{ 
// 		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
// 	}
// 	, 2);

// 	///
// 	/// The actual fitting
// 	///

// 	//Linear fitting
// 	// fitter.SetAlphaWeightingFunction(LinearandFixedBWF);
// 	// fitter.SetBetaWeightingFunction(LinearandFixedBWF);
// 	// double linearandBetaFixedInitialGuess [] = {0.9, 0.2, 0.09, 0.01};
// 	// double* LinearandBetaFixedParams = fitter.Fit(linearandBetaFixedInitialGuess,true);
// 	// //Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kBlue+2);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Linear", H460FittingParamsNewAndImproved, LinearandFixedBWF, LinearandFixedBWF, LinearandBetaFixedParams);

// 	// // //Quadratic fitting
// 	// fitter.SetAlphaWeightingFunction(QuadraticLinearFixedBWF);
// 	// fitter.SetBetaWeightingFunction(LinearandFixedBWF);
// 	// double quadraticAndBetaFixedInitialGuess [] = {0.2, 0., 0.01, 0.12, 0.01};
// 	// double* quadraticAndBetaFixedParams = fitter.Fit(quadraticAndBetaFixedInitialGuess,true);
// 	// // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kRed+2);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, LinearandFixedBWF, quadraticAndBetaFixedParams);

// 	// // //Cubic Fitting
// 	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearFixedBWF);
// 	// fitter.SetBetaWeightingFunction(LinearandFixedBWF);	
// 	// double cubicQuadraticAndBetaFixedInitialGuess [] = {0.06, 0.1, -0.063, 0.002, 0.13, 0.01};
// 	// double* cubicQuadraticAndBetaFixedParams = fitter.Fit(cubicQuadraticAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kGreen+4);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic", H460FittingParamsNewAndImproved, CubicQuadraticLinearFixedBWF, LinearandFixedBWF, cubicQuadraticAndBetaFixedParams);

// 	// //Cubic, no linear offset Fitting
// 	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearBWF);
// 	// fitter.SetBetaWeightingFunction(FixedBWF);	
// 	// double cubicQuadraticLinearAndBetaFixedInitialGuess [] = {0.1177, -0.02, 0.003, 0.13};
// 	// double* cubicQuadraticLinearAndBetaFixedParams = fitter.Fit(cubicQuadraticLinearAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kBlue+4);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic, sans offset", H460FittingParamsNewAndImproved, CubicQuadraticLinearBWF, FixedBWF, cubicQuadraticLinearAndBetaFixedParams);

// 	// //Gaussian Fitting
// 	fitter.SetAlphaWeightingFunction(GaussianBWF);
// 	fitter.SetBetaWeightingFunction(LinearandFixedBWF);	
// 	double gaussianInitialGuess [] = {0.0008, 95.4568, 0.2061, 0.01}; //5, 20, 0.13
// 	double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);


// 	// //Save
// 	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_Fits_Linear_and_quadratic_H460.jpg";
// 	// c->SaveAs((TString)outputName); 


// 	// //Plotting the BWFs themselves
// 	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	// TCanvas* c2 = new TCanvas("c","c");
// 	// c2->SetCanvasSize(4000, 2400);
// 	// // auto legend2 = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
// 	// auto legend2 = new TLegend(0.72,0.14,0.72+0.2,0.14+0.16);//x start, y start, x end, yend
// 	// legend2->SetTextSize(0.04);

// 	// //Setup the marker attributes
// 	// TAttLine lineStyle2{};
	
// 	// lineStyle2.SetLineWidth(5);
// 	// lineStyle2.SetLineStyle(1);

// 	// // //Plot the BWFs
// 	// // lineStyle2.SetLineColor(kGreen+4);
// 	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic", CubicQuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams, "AL", 0., 100.);
// 	// // // lineStyle2.SetLineColor(kBlue+4);
// 	// // // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic, sans offset", CubicQuadraticLinearBWF, cubicQuadraticLinearAndBetaFixedParams, "L", 0., 100.);
// 	// lineStyle2.SetLineColor(kBlue+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Linear", LinearandFixedBWF, LinearandBetaFixedParams,"AL");
// 	// lineStyle2.SetLineColor(kRed+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Quadratic", QuadraticLinearFixedBWF, quadraticAndBetaFixedParams, "L");

// 	// //Save
// 	// std::string outputName2 = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_BWFs_Linear_and_Quadratic_H460.jpg";
// 	// c2->SaveAs((TString)outputName2); 
// }

// void QuadraticBetaFitting()
// {
// 	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
// 	SetupH460SurvivalParameters();

// 	//Plot
// 	//Setup the canvas
// 	gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	TCanvas* c = new TCanvas("c","c");
// 	c->SetCanvasSize(9000, 5000);
// 	c->SetFillStyle(4000);
// 	c->SetFrameFillStyle(4000);
// 	c->Divide(4,3,0.000000005,0.001);
// 	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
// 	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
// 	legend->SetTextSize(0.05);

// 	//Setup the marker attributes
// 	TAttLine lineStyle{};
// 	lineStyle.SetLineColor(kGreen+2);
// 	lineStyle.SetLineWidth(3);
// 	lineStyle.SetLineStyle(2);

// 	//Plot the survival data
// 	SurvivalDataMultigraph(c, legend, H460Params);

// 	//
// 	//Set up the BWFs
// 	//

// 	//Set up the fitter
// 	BWF_Fitter_AlphaBeta fitter{};
// 	fitter.SetCellStudyParameters(H460FittingParamsNewAndImproved);

// 	//Create fixed "BWF"
// 	BiologicalWeightingFunction FixedBWF;
// 	FixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[0];}, 1);

// 	//Create linear BWF
// 	BiologicalWeightingFunction LinearandFixedBWF;
// 	LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

// 	//Exponential growth BWF
// 	BiologicalWeightingFunction ExponentialGrowthBWF;
// 	ExponentialGrowthBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

// 	//Decaying to zero BWF
// 	BiologicalWeightingFunction ExponentialToZeroBWF;
// 	ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

// 	BiologicalWeightingFunction ExponentialToZeroRestrictedBWF;
// 	ExponentialToZeroRestrictedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return 0.15*(1-std::exp(-params[0]/linealEnergy)) ;}, 1);

// 	// //Decaying to zero BWF, simple
// 	// BiologicalWeightingFunction ExponentialToZeroBWF;
// 	// ExponentialToZeroBWF.SetWeightingFunction([](double* params, double linealEnergy) {return params[1]*(1-std::exp(params[0]/linealEnergy)) ;}, 2);

// 	//Simple exponential decay BWF
// 	BiologicalWeightingFunction ExponentialDecayBWF;
// 	ExponentialDecayBWF.SetWeightingFunction([](double* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

// 	//Create quadratic BWF
// 	BiologicalWeightingFunction QuadraticLinearFixedBWF;
// 	QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearFixedBWF;
// 	CubicQuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

// 	//Create cubic BWF
// 	BiologicalWeightingFunction CubicQuadraticLinearBWF;
// 	CubicQuadraticLinearBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy*linealEnergy)+(params[1]*linealEnergy*linealEnergy)+(params[0]*linealEnergy));}, 3);

// 	//Create factored quadratic BWF
// 	BiologicalWeightingFunction QuadraticFactoredBWF;
// 	QuadraticFactoredBWF.SetWeightingFunction([](double* params, double linealEnergy) {return ((params[1]*(linealEnergy-params[2])*(linealEnergy-params[2]))+params[0]);}, 3);

// 	//BiologicalWeightingFunction GaussianBWF;
// 	BiologicalWeightingFunction GaussianBWF;
// 	GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
// 	{ 
// 		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
// 	}
// 	, 2);

// 	///
// 	/// The actual fitting
// 	///

// 	// Linear fitting
// 	fitter.SetAlphaWeightingFunction(LinearandFixedBWF);
// 	fitter.SetBetaWeightingFunction(QuadraticLinearFixedBWF);
// 	double linearandBetaFixedInitialGuess [] = {0.9, 0.2, 0.09, 0.01, 0.001};
// 	double* LinearandBetaFixedParams = fitter.Fit(linearandBetaFixedInitialGuess,true);
// 	//Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kBlue+2);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Linear", H460FittingParamsNewAndImproved, LinearandFixedBWF, QuadraticLinearFixedBWF, LinearandBetaFixedParams);

// 	//Quadratic fitting
// 	fitter.SetAlphaWeightingFunction(QuadraticLinearFixedBWF);
// 	fitter.SetBetaWeightingFunction(QuadraticLinearFixedBWF);
// 	double quadraticAndBetaFixedInitialGuess [] = {0.2, 0., 0.01, 0.12, 0.01, 0.001};
// 	double* quadraticAndBetaFixedParams = fitter.Fit(quadraticAndBetaFixedInitialGuess,true);
// 	// Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kRed+2);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Quadratic", H460FittingParamsNewAndImproved, QuadraticLinearFixedBWF, QuadraticLinearFixedBWF, quadraticAndBetaFixedParams);

// 	//Cubic Fitting
// 	fitter.SetAlphaWeightingFunction(CubicQuadraticLinearFixedBWF);
// 	fitter.SetBetaWeightingFunction(QuadraticLinearFixedBWF);	
// 	double cubicQuadraticAndBetaFixedInitialGuess [] = {0.06, 0.1, -0.063, 0.002, 0.13, 0.01, 0.001};
// 	double* cubicQuadraticAndBetaFixedParams = fitter.Fit(cubicQuadraticAndBetaFixedInitialGuess,true);
// 	// // Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kGreen+4);
// 	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic", H460FittingParamsNewAndImproved, CubicQuadraticLinearFixedBWF, QuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams);

// 	// //Cubic, no linear offset Fitting
// 	// fitter.SetAlphaWeightingFunction(CubicQuadraticLinearBWF);
// 	// fitter.SetBetaWeightingFunction(FixedBWF);	
// 	// double cubicQuadraticLinearAndBetaFixedInitialGuess [] = {0.1177, -0.02, 0.003, 0.13};
// 	// double* cubicQuadraticLinearAndBetaFixedParams = fitter.Fit(cubicQuadraticLinearAndBetaFixedInitialGuess,true);
// 	// // // Plot the BWF prediction on the survival data
// 	// lineStyle.SetLineColor(kBlue+4);
// 	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Cubic, sans offset", H460FittingParamsNewAndImproved, CubicQuadraticLinearBWF, FixedBWF, cubicQuadraticLinearAndBetaFixedParams);

// 	// //Gaussian Fitting
// 	fitter.SetAlphaWeightingFunction(GaussianBWF);
// 	fitter.SetBetaWeightingFunction(QuadraticLinearFixedBWF);	
// 	double gaussianInitialGuess [] = {0.0008, 95.4568, 0.2061, 0.01, 0.001}; //5, 20, 0.13
// 	double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);


// 	// //Save
// 	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_Fits_Linear_and_quadratic_H460.jpg";
// 	// c->SaveAs((TString)outputName); 


// 	// //Plotting the BWFs themselves
// 	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	// TCanvas* c2 = new TCanvas("c","c");
// 	// c2->SetCanvasSize(4000, 2400);
// 	// // auto legend2 = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
// 	// auto legend2 = new TLegend(0.72,0.14,0.72+0.2,0.14+0.16);//x start, y start, x end, yend
// 	// legend2->SetTextSize(0.04);

// 	// //Setup the marker attributes
// 	// TAttLine lineStyle2{};
	
// 	// lineStyle2.SetLineWidth(5);
// 	// lineStyle2.SetLineStyle(1);

// 	// // //Plot the BWFs
// 	// // lineStyle2.SetLineColor(kGreen+4);
// 	// // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic", CubicQuadraticLinearFixedBWF, cubicQuadraticAndBetaFixedParams, "AL", 0., 100.);
// 	// // // lineStyle2.SetLineColor(kBlue+4);
// 	// // // BWFFunctionPlotter(c2, legend2, lineStyle2, "Cubic, sans offset", CubicQuadraticLinearBWF, cubicQuadraticLinearAndBetaFixedParams, "L", 0., 100.);
// 	// lineStyle2.SetLineColor(kBlue+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Linear", LinearandFixedBWF, LinearandBetaFixedParams,"AL");
// 	// lineStyle2.SetLineColor(kRed+2);
// 	// BWFFunctionPlotter(c2, legend2, lineStyle2, "Quadratic", QuadraticLinearFixedBWF, quadraticAndBetaFixedParams, "L");

// 	// //Save
// 	// std::string outputName2 = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Fixed_Beta_BWFs_Linear_and_Quadratic_H460.jpg";
// 	// c2->SaveAs((TString)outputName2); 
// }

// void LQFitting()
// {
// 	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
// 	SetupH460SurvivalParameters();

// 	//Plot
// 	//Setup the canvas
// 	gStyle->SetOptStat(0); //Don't print the stats window in the top right
// 	TCanvas* c = new TCanvas("c","c");
// 	c->SetCanvasSize(9000, 5000);
// 	c->SetFillStyle(4000);
// 	c->SetFrameFillStyle(4000);
// 	c->Divide(4,3,0.000000005,0.001);
// 	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
// 	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
// 	legend->SetTextSize(0.05);

// 	//Setup the marker attributes
// 	TAttLine lineStyle{};
// 	lineStyle.SetLineColor(kGreen+2);
// 	lineStyle.SetLineWidth(3);
// 	lineStyle.SetLineStyle(0);

// 	//Plot the survival data
// 	SurvivalDataMultigraph(c, legend, H460Params);

// 	//Set up the fitter
// 	AlphaBeta_Fitter fitter{};
// 	fitter.SetCellStudyParameters(H460FittingParamsNewAndImproved);

// 	double* AlphaBeta = fitter.Fit(nullptr,true);

// 	//Plot the BWF prediction on the survival data
// 	lineStyle.SetLineColor(kRed+2);
// 	MultigraphSurvivalFitPlotter(c, legend, lineStyle, "", H460FittingParamsNewAndImproved, AlphaBeta);
// 	// SurvivalDataMultigraphResiduals(c, legend, H460FittingParamsNewAndImproved, AlphaBeta);

// 	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/weighted_LQ_fit_H460.jpg";
// 	c->SaveAs((TString)outputName); 
// }

	
void H460Fitting() 
{
	FixedBetaFitting();
	// 	LETFixedBetaFitting();
}
