#pragma once

class TH1F;
//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>

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
#include "TLegend.h"

namespace ErikSpectra
{
	std::vector<std::pair<std::string,TH1D>> GetMonoenergeticLinealSpectra();
	void  PlotMonoenergeticLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& monoenergeticSpectra);
	void  PrintLESSpectra(const std::vector<std::pair<std::string,TH1D>>& monoenergeticSpectra);
};
