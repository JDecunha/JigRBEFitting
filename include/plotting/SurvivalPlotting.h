#pragma once

#include <vector>
class TCanvas; 
class TLegend;
class TH1D;
class CellStudySurvivalParameters;
class CellStudyBWFFittingParameters;
class BWF_Fitting_Results;


// void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudySurvivalParameters survivalParams);
void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudyBWFFittingParameters survivalParams);
void PlotSurvivalFromAlphaBeta(TCanvas* c, const std::vector<std::vector<double>> doselist, const std::vector<double> alphas, const std::vector<double> betas);

void PlotAlphaBeta(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetas, bool plotBeta);
void PlotAlphaBetaFromBWF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, BWF_Fitting_Results results, bool plotBeta);

void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, double const* alphaBetasProton);
void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, BWF_Fitting_Results results);
void PlotRBE10SFLET(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, BWF_Fitting_Results results);
void PlotRBE10SFMcNamara(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium);
void PlotRBE10SFChenAndAhmad(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium);
void PlotRBE10SFWedenberg(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium);