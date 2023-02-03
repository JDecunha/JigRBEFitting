#pragma once

class TH1F;
#include <vector>
#include <string>

namespace ProtonSpectra
{
	//Grabs a single spectrum from specified column
	TH1F GetFadaSpectrumRebinned(const std::string& column);

	//Grabs all the spectra
	std::vector<std::pair<std::string,TH1F>> GetFadaSpectraRebinned();

	void PrintProtonSpectrumForGeantValidation(const std::string& column);
	
	//Saves the spectra to a file
	void SaveFadaSpectra();
};
