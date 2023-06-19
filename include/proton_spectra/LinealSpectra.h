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

namespace LinealSpectra
{
	//Related to the monoenergetic lineal energy spectra
	
		//Get the monoenergetic lineal energy spectra from 0.1-79 MeV (these are the energies in the KE spectrum from Fada)
		std::vector<std::pair<std::string,TH1D>> GetMonoenergeticLinealSpectra();
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotMonoenergeticLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& monoenergeticSpectra);
		
	//Related to the KE weighted lineal energy spectra. For the different locations in the jig.

		//Get the lineal energy for each of Fada's cell samples
		std::vector<std::pair<std::string,TH1D>> GetKeWeightedFrequencyLinealSpectra(std::string targetSize);
		//This is for dose spectra
		std::vector<std::pair<std::string,TH1D>> GetKeWeightedLinealSpectra(std::string targetSize);
		//New proton spectrum
		std::vector<std::pair<std::string,TH1D>> GetKeWeightedLinealSpectraJuly2023(std::string targetSize);
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotKeWeightedLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& linealLibrary);
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotKeWeightedLinealSpectraMultigraph(const std::vector<std::pair<std::string,TH1D>>& linealLibrary);
		
		//Save the weighted lineal spectra
		void SaveKeWeightedLinealSpectra();
		void SaveKeWeightedFrequencyLinealSpectra();
		void SaveKeWeightedLinealSpectraCSV();

		//For the Lineal energy library manuscript
		void CompareGeant4toFluenceWeighting();
};
