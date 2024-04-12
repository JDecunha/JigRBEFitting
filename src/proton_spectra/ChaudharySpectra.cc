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
#include "ChaudharySpectra.h"
#include "Utilities.h"

namespace utils = Utilities;

std::vector<std::pair<std::string,TH1F>> ChaudharySpectra::GetChaudharyKESpectra(std::string filePath, std::vector<std::string> columns)
{
	std::vector<std::pair<std::string,TH1F>> output;
	TH1::AddDirectory(kFALSE); //so that the TH1s don't get added to the global gDirectory by default
	TTree KEspectra = TTree("T","a collection of KE spectra from Chaudhary et al. 2014"); //Laod the file into a TTree
	KEspectra.ReadFile((TString)filePath); //all variables will be floating point

	for (const auto& column:columns)
	{
		//Make the TTreeReader to read the file
		TTreeReader treeReader(&KEspectra);
		TTreeReaderValue<Float_t> KeBin(treeReader,"Energy"); //Energy bins in the leftmost column
		TTreeReaderValue<Float_t> KeValue(treeReader,(TString)column); //This will pull out the column for the LET currently selected

		//Make a vector to hold the bins
		std::vector<double> KE_bins;

		//Put the energy bins into the std::vector
		while (treeReader.Next()) 
		{
			KE_bins.push_back(*KeBin);
			// std::cout << *KeBin << std::endl;
		}

		//Use this information to fill a TH1F
		TH1F KE_spectrum = TH1F("KE Spectrum"," ",(KE_bins.size()-1), KE_bins.data());
		treeReader.Restart(); //Restart the reader so it can run again again

		//Fill the TH1F with values
		while (treeReader.Next()) 
		{
			KE_spectrum.Fill(*KeBin,*KeValue); 
		}

		//For the Chaudhary data we Don't need to rebin because I did it already in Libreoffice Calc

		//Part 2: We shouldn't need to normalize the spectrum either, because LibreOffice Calc says they're normalized
		//But we will anyways
		float normalization = 0;

		for (int i = 0; i <= KE_spectrum.GetNbinsX(); ++i) //calculate the normalization value
		{
			normalization += KE_spectrum.GetBinContent(i);
		}

		for (int i = 0; i <= KE_spectrum.GetNbinsX(); ++i) //normalize
		{
			float normalizedVal = KE_spectrum.GetBinContent(i)/normalization;
			KE_spectrum.SetBinContent(i,normalizedVal);
			// std::cout << KE_spectrum.GetBinCenter(i) << " " << KE_spectrum.GetBinContent(i) << std::endl;
		}

		output.push_back(std::pair(column, KE_spectrum));

	}

	return output;
}


void SaveChaudharyKESpectra()
{
	auto chaudharyPristine = ChaudharySpectra::GetChaudharyKESpectra("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/Pristine_rebin.csv", std::vector<std::string> {"1.1","3.9","6.7","11.6","17.7","22.5"});
	auto chaudharySOBP = ChaudharySpectra::GetChaudharyKESpectra("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/SOBP_rebin.csv", std::vector<std::string> {"1.27","3","4.4","13.7","20.9","25.4"});

	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/ChaudharyPristineKESpectra.root","RECREATE");
	spectrumFile->WriteObject(&chaudharyPristine, "Spectrum_Library");
	delete spectrumFile;

	spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/ChaudharySOBPKESpectra.root","RECREATE");
	spectrumFile->WriteObject(&chaudharySOBP, "Spectrum_Library");
	delete spectrumFile;
}

std::vector<std::pair<std::string,TH1D>> ChaudharySpectra::GetFy(std::string filePath, std::string targetSize)
{
	/*
	Retrieves the d(y) spectra calculated in the chaudhary library
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/SuperTrack_01to100MeV_June2023_binspanning";

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<std::string,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open((TString)filePath);
	spectrumFile->GetObject("Spectrum_Library",keSpectra);
	spectrumFile->Close();
	
	//This will hold the output
	std::vector<std::pair<std::string,TH1D>> linealEnergyLibrary; 
	
	for(const auto& spectra:*keSpectra)
	{
		//Step 1: Weight the lineal energy spectrum by the KE spectrum.
		//Grab the KE specta for a given irradiation
		TH1F KESpectrum = std::get<1>(spectra);
		std::string KEColumnName = std::get<0>(spectra);

		//Just get the bins for the output spectrum
		TH1D outputLinealSpectrum = utils::GetDy(fyFolder, 0.15, targetSize);
		outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

		//Loop over the KE spectrum and weight the d(y) spectra accordingly
		for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			TH1D tempLinealSpectrum = utils::GetNy(fyFolder, KESpectrum.GetBinCenter(i), targetSize); //Get the N(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Scale the spectrum so N(y) matches an effective count of 10 million
			long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, KESpectrum.GetBinCenter(i), targetSize);
			double scalingFactor = double(1e7)/double(effectiveNumTracks);
			tempLinealSpectrum.Scale(scalingFactor);

			//Add the weighted spectrum to the output fy spectrum
			outputLinealSpectrum.Add(&tempLinealSpectrum);
		}	
		
		//Convert N(y) to d(y)
		utils::PMF_to_FrequencyFunction(&outputLinealSpectrum);
		//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7
		linealEnergyLibrary.push_back(std::make_pair<std::string,TH1D>(std::move(KEColumnName),std::move(outputLinealSpectrum)));
	}

	delete keSpectra; //we own the keSpectra pointer so we have to delete it

	return linealEnergyLibrary;
}

std::vector<std::pair<std::string,TH1D>> ChaudharySpectra::GetDy(std::string filePath, std::string targetSize)
{
	/*
	Retrieves the d(y) spectra calculated in the chaudhary library
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/SuperTrack_01to100MeV_June2023_binspanning";

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<std::string,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open((TString)filePath);
	spectrumFile->GetObject("Spectrum_Library",keSpectra);
	spectrumFile->Close();
	
	//This will hold the output
	std::vector<std::pair<std::string,TH1D>> linealEnergyLibrary; 
	
	for(const auto& spectra:*keSpectra)
	{
		//Step 1: Weight the lineal energy spectrum by the KE spectrum.
		//Grab the KE specta for a given irradiation
		TH1F KESpectrum = std::get<1>(spectra);
		std::string KEColumnName = std::get<0>(spectra);

		//Just get the bins for the output spectrum
		TH1D outputLinealSpectrum = utils::GetDy(fyFolder, 0.15, targetSize);
		outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

		//Loop over the KE spectrum and weight the d(y) spectra accordingly
		for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			TH1D tempLinealSpectrum = utils::GetNy(fyFolder, KESpectrum.GetBinCenter(i), targetSize); //Get the N(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Scale the spectrum so N(y) matches an effective count of 10 million
			long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, KESpectrum.GetBinCenter(i), targetSize);
			double scalingFactor = double(1e7)/double(effectiveNumTracks);
			tempLinealSpectrum.Scale(scalingFactor);

			//Add the weighted spectrum to the output fy spectrum
			outputLinealSpectrum.Add(&tempLinealSpectrum);
		}	
		
		//Convert N(y) to d(y)
		utils::PMF_to_DoseFunction(&outputLinealSpectrum);
		//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7
		linealEnergyLibrary.push_back(std::make_pair<std::string,TH1D>(std::move(KEColumnName),std::move(outputLinealSpectrum)));
	}

	delete keSpectra; //we own the keSpectra pointer so we have to delete it

	return linealEnergyLibrary;
}

void SaveFySpectraCSV(std::string filePath, std::string outputPath, std::string outputName)
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& targetSize : targetSizes)
	{
		auto keWeightedSpectra = ChaudharySpectra::GetFy(filePath,targetSize);

		for (const auto& linealSpectra:keWeightedSpectra)
		{
			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			auto columnName = std::get<0>(linealSpectra);

			std::ofstream myfile;
			myfile.open (outputPath+"/"+outputName+columnName+".csv");
			int length = spectra->GetNbinsX(); //get length

			for (int i = 1; i <= length; i++) 
			{
				auto value = spectra->GetBinContent(i);
				double low_edge = spectra->GetBinLowEdge(i);
				myfile << std::setprecision(6) << low_edge << "," <<  value <<"\n";

			}

			myfile.close();
		}

	}
}

void SaveDySpectraCSV(std::string filePath, std::string outputPath, std::string outputName)
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& targetSize : targetSizes)
	{
		auto keWeightedSpectra = ChaudharySpectra::GetDy(filePath,targetSize);

		for (const auto& linealSpectra:keWeightedSpectra)
		{
			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			auto columnName = std::get<0>(linealSpectra);

			std::ofstream myfile;
			myfile.open (outputPath+"/"+outputName+columnName+".csv");
			int length = spectra->GetNbinsX(); //get length

			for (int i = 1; i <= length; i++) 
			{
				auto value = spectra->GetBinContent(i);
				double low_edge = spectra->GetBinLowEdge(i);
				myfile << std::setprecision(6) << low_edge << "," <<  value <<"\n";

			}

			myfile.close();
		}

	}
}

void ChaudharySaveLinealSpectraWrapper()
{
	SaveFySpectraCSV("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/ChaudharyPristineKESpectra.root","/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/fy","PristineF1y");
	SaveFySpectraCSV("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/ChaudharySOBPKESpectra.root","/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/fy","SOBPFy");

	SaveDySpectraCSV("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/ChaudharyPristineKESpectra.root","/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/dy","PristineDy");
	SaveDySpectraCSV("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/ChaudharySOBPKESpectra.root","/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/Chaudhary/dy","SOBPDy");
}


void ChaudharySpectra::main()
{
	std::cout << "we straight chillin" << std::endl;
	ChaudharySaveLinealSpectraWrapper();
}