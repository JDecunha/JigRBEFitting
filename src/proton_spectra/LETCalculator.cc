//ROOT
#include "TH1F.h"
//std
#include <fstream>
#include <sstream>
#include <string>
//This project
#include "ProtonSpectra.h"
#include "LinealSpectra.h"
#include "Utilities.h"
#include "LETCalculator.h"

std::pair<std::vector<double>,std::vector<double>> PullStoppingPowersFromCsv()
{
	//The two vectors that get filled
	std::vector<double> protonEnergy;
	std::vector<double> protonElecStoppingPower;

	//The path, an input stream, and a temp to hold the line being readout
	std::string filePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Documents/NIST_CSV.csv";
	std::ifstream SPFile(filePath);
	std::string ReadoutLine;

	while (std::getline(SPFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		
		//Get the two comma separated values and fill
		std::getline(spaceSplitter,word, ',');
		protonEnergy.push_back(std::stod(word));
		std::getline(spaceSplitter,word, ',');
		protonElecStoppingPower.push_back(std::stof(word));
	}

	return std::make_pair(protonEnergy,protonElecStoppingPower);
}


std::pair<std::vector<double>,std::vector<double>> PullStoppingPowersFromCsvAndInterpolate()
{
	//Step 1.) Fill data from file

	//The two vectors that get filled
	std::vector<double> protonEnergy;
	std::vector<double> protonElecStoppingPower;

	//The path, an input stream, and a temp to hold the line being readout
	std::string filePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Documents/NIST_PSTAR_FULL.csv";
	std::ifstream SPFile(filePath);
	std::string ReadoutLine;

	while (std::getline(SPFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		
		//Get the two comma separated values and fill
		std::getline(spaceSplitter,word, ',');
		protonEnergy.push_back(std::stod(word));
		std::getline(spaceSplitter,word, ',');
		protonElecStoppingPower.push_back(std::stof(word));
	}

	//Step 2.) Interpolate with 0.1 MeV spacing
	std::vector<double> interpolatedPstarEnergies = Utilities::Linspace(0.1,80,800); //These are all the bin values we want
	std::vector<double> protonElecStoppingPowerInterp;

	// for (auto val:interpolatedPstarEnergies)
	// {
	// 	std::cout << val << std::endl;
	// }

	for (int i = 0; i < interpolatedPstarEnergies.size(); ++i) //Loop over the desired interpolation values
	{
		double lowPstarBin = 0.;
		double highPstarBin = 0.; 

		for (int j = 0; j < protonEnergy.size(); ++j) //Then loop over the regular PStar values
		{
			if (interpolatedPstarEnergies[i] == protonEnergy[j]) //If the interpolation energy matches the PStar energy, just take it
			{
				protonElecStoppingPowerInterp.push_back(protonElecStoppingPower[j]);
				break;
			}

			if (protonEnergy[j] > interpolatedPstarEnergies[i]) //If we have passed the bin that contains the energy we want
			{
				lowPstarBin = protonEnergy[j-1]; //We passed the desired lower bin energy, so grab the last one from pstar
				highPstarBin = protonEnergy[j];

				double slope = (protonElecStoppingPower[j]-protonElecStoppingPower[j-1])/(highPstarBin-lowPstarBin);
				double distance = interpolatedPstarEnergies[i]-lowPstarBin;

				double interpolatedStoppingPower = (protonElecStoppingPower[j-1])+(slope*distance);

				protonElecStoppingPowerInterp.push_back(interpolatedStoppingPower);
				break;
			}

		}

	}

	return std::make_pair(interpolatedPstarEnergies,protonElecStoppingPowerInterp);
}


void ComputeLETs()
{
	// Step 1.) Get the Stopping Power Values
	std::pair<std::vector<double>,std::vector<double>> NISTPSTAR = PullStoppingPowersFromCsv();

	// Step 2.) Pull all of the fluence spectra
	std::vector<std::pair<std::string,TH1F>> spectra = ProtonSpectra::GetFadaSpectraRebinned();	

	// Step 3.) Compute LET for each fluence spectrum
	for (auto& pair:spectra) //For every set of fluence spectra that we have
	{

		TH1F spectrum = std::get<1>(pair); //Grab the fluence spectrum
		double LETTopSummation = 0; double LETBottomSummation = 0; //Accumulators

		for (int i = 2; i <= spectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{

			//Because of the way we re-defined the fluence spectrum, "GetBinLowEdge" actually corresponds to the middle
			double keSpectrumEnergy = spectrum.GetBinLowEdge(i); double keSpectrumBinValue = spectrum.GetBinContent(i);

			for (int pStarLoopIndex = 0; pStarLoopIndex < std::get<0>(NISTPSTAR).size(); ++pStarLoopIndex) //Loop over PSTAR to find the Stopping power
			{

				double pStarEnergyValue = std::get<0>(NISTPSTAR)[pStarLoopIndex];

				if (std::fabs(pStarEnergyValue - keSpectrumEnergy) < 0.001) //We found entry in PSTAR if true
				{

					double pStarStoppingPower = std::get<1>(NISTPSTAR)[pStarLoopIndex];

					//Add to the accumulators
					LETTopSummation += pStarStoppingPower*pStarStoppingPower*keSpectrumBinValue;
					LETBottomSummation += pStarStoppingPower*keSpectrumBinValue;

					// std::cout << "Match found" << keSpectrumEnergy << ", " << pStarEnergyValue << std::endl;

					break; //We found it, we can end the for loop

				}
			}
		}	

		//Accumulators are done, actually sum LET
		double LETd = LETTopSummation*0.1/LETBottomSummation;
		std::cout << std::get<0>(pair) << ", LETd: " << LETd << std::endl;	
	}

}

void ComputeLETsTakeTwo()
{
	// Step 1.) Get the Stopping Power Values
	std::pair<std::vector<double>,std::vector<double>> NISTPSTAR = PullStoppingPowersFromCsv();

	// Step 2.) Pull all of the fluence spectra
	std::vector<std::pair<std::string,TH1F>> spectra = ProtonSpectra::GetFadaKESpectrumRaw();	

	// Step 3.) Compute LET for each fluence spectrum
	for (auto& pair:spectra) //For every set of fluence spectra that we have
	{

		TH1F spectrum = std::get<1>(pair); //Grab the fluence spectrum
		double LETTopSummation = 0; double LETBottomSummation = 0; //Accumulators

		for (int i = 1; i < spectrum.GetNbinsX(); ++i) //start at 1 to skip underflow
		{

			//Because we haven't rebinnned Fada's spectrum, we interpolate the bins to get the midpoints
			double keSpectrumEnergy = (spectrum.GetBinCenter(i)+spectrum.GetBinCenter(i+1))/2.; 
			// std::cout << keSpectrumEnergy << std::endl;
			double keSpectrumBinValue = (spectrum.GetBinContent(i)+spectrum.GetBinContent(i+1))/2.;

			for (int pStarLoopIndex = 0; pStarLoopIndex < std::get<0>(NISTPSTAR).size(); ++pStarLoopIndex) //Loop over PSTAR to find the Stopping power
			{

				double pStarEnergyValue = std::get<0>(NISTPSTAR)[pStarLoopIndex];

				if (std::fabs(pStarEnergyValue - keSpectrumEnergy) < 0.001) //We found entry in PSTAR if true
				{

					double pStarStoppingPower = std::get<1>(NISTPSTAR)[pStarLoopIndex];

					//Add to the accumulators
					LETTopSummation += pStarStoppingPower*pStarStoppingPower*keSpectrumBinValue;
					LETBottomSummation += pStarStoppingPower*keSpectrumBinValue;

					// std::cout << "Match found" << keSpectrumEnergy << ", " << pStarEnergyValue << std::endl;

					break; //We found it, we can end the for loop

				}
			}
		}	

		//Accumulators are done, actually sum LET
		double LETd = LETTopSummation*0.1/LETBottomSummation;
		std::cout << std::get<0>(pair) << ", LETd: " << LETd << std::endl;	
	}

};

void ComputeLETsTakeThree()
{
	// Step 1.) Get the Stopping Power Values
	std::pair<std::vector<double>,std::vector<double>> NISTPSTAR = PullStoppingPowersFromCsvAndInterpolate();

	// Step 2.) Pull all of the fluence spectra
	std::vector<std::pair<std::string,TH1F>> spectra = ProtonSpectra::GetFadaKESpectrumRaw();	

	// Step 3.) Compute LET for each fluence spectrum
	for (auto& pair:spectra) //For every set of fluence spectra that we have
	{

		TH1F spectrum = std::get<1>(pair); //Grab the fluence spectrum
		double LETTopSummation = 0; double LETBottomSummation = 0; //Accumulators

		for (int i = 1; i < spectrum.GetNbinsX(); ++i) //start at 1 to skip underflow
		{

			//Because we haven't rebinnned Fada's spectrum, we interpolate the bins to get the midpoints
			double keSpectrumEnergy = (spectrum.GetBinCenter(i)+spectrum.GetBinCenter(i+1))/2.; 
			// std::cout << keSpectrumEnergy << std::endl;
			double keSpectrumBinValue = (spectrum.GetBinContent(i)+spectrum.GetBinContent(i+1))/2.;

			for (int pStarLoopIndex = 0; pStarLoopIndex < std::get<0>(NISTPSTAR).size(); ++pStarLoopIndex) //Loop over PSTAR to find the Stopping power
			{

				double pStarEnergyValue = std::get<0>(NISTPSTAR)[pStarLoopIndex];

				if (std::fabs(pStarEnergyValue - keSpectrumEnergy) < 0.001) //We found entry in PSTAR if true
				{

					double pStarStoppingPower = std::get<1>(NISTPSTAR)[pStarLoopIndex];

					//Add to the accumulators
					LETTopSummation += pStarStoppingPower*pStarStoppingPower*keSpectrumBinValue;
					LETBottomSummation += pStarStoppingPower*keSpectrumBinValue;

					// std::cout << "Match found" << keSpectrumEnergy << ", " << pStarEnergyValue << std::endl;

					break; //We found it, we can end the for loop

				}
			}
		}	

		//Accumulators are done, actually sum LET
		double LETd = LETTopSummation*0.1/LETBottomSummation;
		std::cout << std::get<0>(pair) << ", LETd: " << LETd << std::endl;	
	}

};

void LETCalcMain()
{
	// std::cout << "Take  one: " << std::endl;
	// ComputeLETs();
	// std::cout << "\nTake two: " << std::endl;
	// ComputeLETsTakeTwo();

	ComputeLETsTakeThree();
}