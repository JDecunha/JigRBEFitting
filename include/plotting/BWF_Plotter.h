#pragma once

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
#include "SurvivalPlotting.h"
#include "BWF_Fitting_Results.h"
#include "NonLinearLeastSquaresUtilities.h"

//Functions for fitting and plotting
TGraph* BWFFunctionPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, BiologicalWeightingFunction fittingFunction, double* fitFuncParams, std::string options, double minLineal = 0, double maxLineal = 250);
void GeneralizedBWFMultigraphPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BWF_Fitting_Results results, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);
void AlphaBetaMultigraphResiduals(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BWF_Fitting_Results results,  double minDose = 0, double maxDose = 5.75);

//
// Depreceated
//
void GeneralizedBWFMultigraphPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudySurvivalParameters& survivalParams, BiologicalWeightingFunction fittingFunction, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);
void GeneralizedBWFMultigraphPlotterBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudySurvivalParameters& survivalParams, BiologicalWeightingFunction fittingFunction, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);
void GeneralizedBWFMultigraphPlotterAlphaBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);
void GeneralizedBWFMultigraphPlotterAlphaBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BWF_Fitting_Results results,  double minDose = 0, double maxDose = 5.75);

void GeneralizedLETMultigraphPlotterAlphaBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);

void AlphaBetaMultigraphResiduals(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose = 0, double maxDose = 5.75);

void MultigraphSurvivalFitPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, double* alphas_and_betas, double minDose = 0, double maxDose = 5.75);
void SurvivalDataMultigraphResiduals(TCanvas* c, TLegend* legend, CellStudyBWFFittingParameters& survivalParams, double* alphas_and_betas);