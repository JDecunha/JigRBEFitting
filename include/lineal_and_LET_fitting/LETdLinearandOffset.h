#pragma once

class TCanvas; 
class TLegend;
class TH1D;

//GNU Scientific Library
#include <gsl/gsl_vector.h>

int LinearSurvivalLETd(const gsl_vector* x, void* data, gsl_vector* f);
TCanvas* PlotSurvivalFromLinearLETModel(TCanvas* c, TLegend* legend, const std::vector<std::vector<double>> doselist, const std::vector<double> LETs, const double c0, const double c1, const double beta);