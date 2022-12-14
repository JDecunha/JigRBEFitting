#pragma once

class TH1F;

namespace LinealSpectra
{
	//Related to the monoenergetic lineal energy spectra
	
		//Get the monoenergetic lineal energy spectra from 0.1-79 MeV (these are the energies in the KE spectrum from Fada)
		std::vector<std::pair<std::string,TH1D>> GetMonoenergeticLinealSpectra();
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotMonoenergeticLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& monoenergeticSpectra);
		
	//Related to the KE weighted lineal energy spectra. For the different locations in the jig.

		//Get the lineal energy for each of Fada's cell samples
		std::vector<std::pair<std::string,TH1D>> GetKeWeightedLinealSpectra();
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotKeWeightedLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& linealLibrary);
		
		//Plot the weighted lineal energy spectra on a multigraph
		void PlotKeWeightedLinealSpectraMultigraph(const std::vector<std::pair<std::string,TH1D>>& linealLibrary);
		
		//Save the weighted lineal spectra
		void SaveKeWeightedLinealSpectra();
};
