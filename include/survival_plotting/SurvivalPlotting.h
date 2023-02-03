#pragma once

#include <vector>
class TCanvas; 
class TLegend;
class TH1D;
class CellStudySurvivalParameters;


void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudySurvivalParameters survivalParams);
void PlotSurvivalFromAlphaBeta(TCanvas* c, const std::vector<std::vector<double>> doselist, const std::vector<double> alphas, const std::vector<double> betas);