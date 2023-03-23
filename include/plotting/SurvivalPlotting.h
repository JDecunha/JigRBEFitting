#pragma once

#include <vector>
class TCanvas; 
class TLegend;
class TH1D;
class CellStudySurvivalParameters;
class CellStudyBWFFittingParameters;
class BWF_Fitting_Results;


void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudySurvivalParameters survivalParams);
void PlotSurvivalFromAlphaBeta(TCanvas* c, const std::vector<std::vector<double>> doselist, const std::vector<double> alphas, const std::vector<double> betas);

void PlotAlphaBeta(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetas, bool plotBeta);
void PlotAlphaBetaFromBWF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, BWF_Fitting_Results results, bool plotBeta);