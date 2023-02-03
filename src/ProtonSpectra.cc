//std
#include <iostream>
#include <fstream>
#include <filesystem>
#include <utility>

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

//This project
#include "ProtonSpectra.h"
#include "Utilities.h"

namespace utils = Utilities;


TH1F ProtonSpectra::GetFadaSpectrumRebinned(const std::string& column)
{
	TH1::AddDirectory(kFALSE); //so that the TH1s don't get added to the global gDirectory by default

	//Part 1.) Read Fada's data in and get it into a TH1F

	//Laod the file into a TTree
	TTree KEspectra = TTree("T","a collection of proton KE spectra");
	KEspectra.ReadFile("/home/joseph/Documents/PHD_local/fada_data/protonKEspectra.csv"); //all variables will be floating point

	//Make the TTreeReader
	TTreeReader treeReader(&KEspectra);
	TTreeReaderValue<Float_t> KeBin(treeReader,"KE_MeV");
	TTreeReaderValue<Float_t> KeValue(treeReader,(TString)column);

	//Make a vector to hold the bins
	std::vector<double> KE_bins;

	//Put the bins into the std::vector
	while (treeReader.Next()) 
	{
		KE_bins.push_back(*KeBin);
	}

	//Make the TH1F with correct bins, it's currently empty
	TH1F KE_spectrum = TH1F("KE Spectrum"," ",(KE_bins.size()-1), KE_bins.data());
	treeReader.Restart(); //Restart the reader so it can run again again

	//Fill the TH1F with values
	while (treeReader.Next()) 
	{
		KE_spectrum.Fill(*KeBin,*KeValue); 

	}

	//Part 2.) Rebinning
	//(Sorry to my future self that I find this very confusing)

	/*
	So, for rebinning, from 0.1-0.9 MeV, we just take the bin above and bin below and average them to get the bin "middle" best guess.
	Since we're not actually creating or combining any new bins we want to keep the normalization as is.

	It's okay that the 0.0 MeV bin only gets half counted, because we can consider half of it to be in the 0 MeV bin which we ignore.

	For 1 MeV, we're taking the average of 0.9 - 1.4 MeV. 
	For this, after the 0.9 MeV bin value in the new spectrum is calculated, we take the 0.9 MeV value from the original spectrum and halve its value. 
	Then when we average from 0.9 MeV - 1.4 MeV we can multiply by 6 and preserve the same number of the total counts. 
	(We are normalizing with the principle in mind, that say 0.9 MeV had 2 counts and the rest of the spectrum zero counts. 
	When we perform averaging there should still be a total of two counts. I.e. we are preserving the total number of counts in the averaging.
	Except for the 0 MeV bin, for reasons explained above.
	You can convince yourself this works by working with the above example of 2 counts in 0.9 MeV and none anywhere else.)

	For the Bins from 2-78 MeV we're taking 10 bins total and combining them into one. So the average should be multiplied by ten to keep the same number of total counts.
	*/

	//Create the newBins we want to output with
	std::vector<double> newBins;
	newBins.reserve(90); //pre-allocate memory so it's faster
	std::vector<double> firstBinValues = utils::Linspace(0,0.9,10); //from 0-0.9 MeV
	std::vector<double> lastBinValues = utils::Linspace(1,80,80); //from 1-80 MeV
	//Insert the new bins into the vector
	newBins.insert(newBins.end(),firstBinValues.begin(),firstBinValues.end()); 
	newBins.insert(newBins.end(),lastBinValues.begin(),lastBinValues.end());

	//New histogram
	TH1F KE_spectrum_averaged = TH1F("KE Spectrum Rebinned"," ",(newBins.size()-1), newBins.data());
	
	//Part 2.1: rebin from 0.1 - 0.9 MeV with averaging
	for (int i = 1; i < 10; i++) 
	{
		std::pair<float,float> binAverageAndCenter = utils::GetHistogramAverage(KE_spectrum, i, i+1);
		KE_spectrum_averaged.Fill(std::get<1>(binAverageAndCenter), std::get<0>(binAverageAndCenter));
	}

	//Part 2.2: rebin 1 MeV. 1 MeV is the 'in between' bin, so average from 0.9-1.4 MeV
	KE_spectrum.SetBinContent(10, KE_spectrum.GetBinContent(10)*0.5); //Halve the 0.9 MeV bin
	std::pair<float,float> binAverageAndCenter = utils::GetHistogramAverage(KE_spectrum, 10, 15); //10th bin is 0.9, 15th is 1.4 MeV
	KE_spectrum_averaged.Fill(std::get<1>(binAverageAndCenter), std::get<0>(binAverageAndCenter)*6); //multiply by the 6 combined bins

	//Part 2.3: rebin from 2-79 MeV
	for (int i = 0; i < 78; ++i) 
	{
		std::pair<float,float> binAverageAndCenter = utils::GetHistogramAverage(KE_spectrum, 16+(10*i), 16+9+(10*i)); //I hate using hardcoded numbers

		if (std::get<1>(binAverageAndCenter) <= 79)
		{
			KE_spectrum_averaged.Fill(i+2, std::get<0>(binAverageAndCenter)*10);
		}
		else {std::cout << "bin out of range" << std::endl;}
	}
 
	//Part 3: Normalizing the spectrum (so it adds to 1)
	float normalization = 0;

	for (int i = 0; i <= KE_spectrum_averaged.GetNbinsX(); ++i) //calculate the normalization value
	{
		normalization += KE_spectrum_averaged.GetBinContent(i);
	}

	for (int i = 0; i <= KE_spectrum_averaged.GetNbinsX(); ++i) //normalize
	{
		float normalizedVal = KE_spectrum_averaged.GetBinContent(i)/normalization;
		KE_spectrum_averaged.SetBinContent(i,normalizedVal);
		//std::cout << KE_spectrum_averaged.GetBinCenter(i) << " " << KE_spectrum_averaged.GetBinContent(i) << std::endl;
	}

	return KE_spectrum_averaged;
}

std::vector<std::pair<std::string,TH1F>> ProtonSpectra::GetFadaSpectraRebinned() //gets all of Fadas spectra and outputs
{
	std::vector<std::pair<std::string,TH1F>> output;
	std::vector<std::string> columns{"a","b","c","d","e","f","g","h","i","j","k","l"};

	for (const auto& item:columns)
	{
		output.push_back(std::pair<std::string,TH1F>(item, GetFadaSpectrumRebinned(item)));
	}

	return output;
}

void ProtonSpectra::SaveFadaSpectra()
{
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");
	auto spectra = ProtonSpectra::GetFadaSpectraRebinned();
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/FadasProtonSpectra.root","RECREATE");
	spectrumFile->WriteObject(&spectra, "Spectrum_Library");
	delete spectrumFile;
}

void ProtonSpectra::PrintProtonSpectrumForGeantValidation(const std::string& column)
{
	TH1F protonSpectrum = GetFadaSpectrumRebinned(column);

	std::filesystem::path outputPath = std::filesystem::path("/home/joseph/Documents/PHD_local/GeantValidationSpectra");
	std::ofstream outputStream(outputPath.string() + "/spectrum_" + column + "_.txt");

	if (outputStream.is_open())
	{
		for (int i = 1; i <= protonSpectrum.GetNbinsX(); ++i)
		{
			outputStream << "/gps/hist/point " << protonSpectrum.GetBinCenter(i) << " " << protonSpectrum.GetBinContent(i) << "\n";
		}
	}
	else {std::cout << "file failed to open" << std::endl;}
	
}

