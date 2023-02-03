#pragma once

//std
#include <vector>
#include <string>
//ROOT
#include "TH1D.h"
//GNU Scientific Library
#include <gsl/gsl_multifit_nlinear.h>
//This project
#include "BiologicalWeightingFunction.h"

class CellStudySurvivalParameters
{
  public:
    //Physical information
    std::vector<std::pair<std::string,TH1D>> dySpectra;
    std::vector<double> LETd;

    //Survival information
    std::vector<std::vector<double>> survivingFraction;
    std::vector<std::vector<double>> survivingFractionUncertainty;
    std::vector<std::vector<double>> dose;

    //Beta from photons
    double beta;

    //The desired fitting function
    BiologicalWeightingFunction fittingFunction;

    int GetNumDataPoints()
    {
      int i = 0;
      for (const auto& experimentlist:dose)
      {
        for (const auto& individualirradiation:experimentlist)
        {
          ++i;
        }
      }
      return i;
    };

};

struct CellStudySurvivalParametersLETd
{
  std::vector<double> LETd;
  std::vector<std::vector<double>> survivingFraction;
  std::vector<std::vector<double>> dose;
  double beta;
};

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);