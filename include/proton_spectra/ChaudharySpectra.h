#pragma once

class TH1F;
#include <vector>
#include <string>

namespace ChaudharySpectra
{
	//Grabs a single spectrum from specified column
	std::vector<std::pair<std::string,TH1F>> GetChaudharyKESpectra(std::string filePath, std::vector<std::string> columns);
	std::vector<std::pair<std::string,TH1D>> GetFy(std::string filePath, std::string targetSize);
	std::vector<std::pair<std::string,TH1D>> GetDy(std::string filePath, std::string targetSize);
	void main();
};
