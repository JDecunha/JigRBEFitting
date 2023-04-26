//std
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <filesystem>
#include <utility>
//CERN ROOT
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
//This project
#include "Uwes_Spectrum.h"
#include "Utilities.h"

namespace utils = Utilities;

//Alright, Uwe's spectrum
//10 keV energy resolution. From 0.01 MeV up to 100 MeV (10 k bins)
//1 mm spatial resolution. From 1 mm to 70 mm (70 spectra)

std::pair<std::vector<double>,std::vector<TH1F>> Import_Spectrum_and_Rebin()
{
	//The path, an input stream, and a temp to hold the line being readout
	std::string filePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Documents/Uwes_First_Spectrum.dat";
	std::ifstream spectrumFile(filePath);
	std::string ReadoutLine;

	//Create doubles corresponding to Uwe's bins
	std::vector<double> uwesBinValues = utils::Linspace(0.01,100, 10000);

	//Create the newBins we want to output with
	//NOTE TO SELF: Keep in mind, this histogram is defined "function style"
	//So even though the first bin "spans" 0-0.1 MeV, it really is representative of the value AT 0.1 MeV
	//i.e. it's not a true histogram
	std::vector<double> newBins;
	newBins.reserve(90); //pre-allocate memory so it's faster
	std::vector<double> firstBinValues = utils::Linspace(0,0.9,10); //from 0-0.9 MeV
	std::vector<double> lastBinValues = utils::Linspace(1,80,80); //from 1-80 MeV
	//Insert the new bins into the vector
	newBins.insert(newBins.end(),firstBinValues.begin(),firstBinValues.end()); 
	newBins.insert(newBins.end(),lastBinValues.begin(),lastBinValues.end());

	//The output of this function
	std::vector<double> distance = utils::Linspace(0.05,70.5, 70);
	std::vector<TH1F> spectra;

	while (std::getline(spectrumFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		int loopValue = 0; int innerLoopValue = 0; 

		//Create a TH1D for the output values
		TH1F Uwes_KE_spectrum_averaged = TH1F("Uwes KE Spectrum Rebinned"," ",(newBins.size()-1), newBins.data());

		//Re-bin from 0 - 0.15 MeV (0.1 MeV bin)
		while (std::getline(spaceSplitter,word, ' '))
		{
			if (loopValue == 0) {std::cout << "first bin val: " << word << std::endl;}
			double valueOver15 = std::stod(word)/15;
			Uwes_KE_spectrum_averaged.Fill(0.05,valueOver15); //0.05 corresponds to 0.1 MeV bin
			// std::cout << "Value: " << word << " Averaged: " << valueOver15 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(1) << std::endl;
			// std::cout << uwesBinValues[loopValue] << std::endl;
			if (std::fabs(uwesBinValues[loopValue] - 0.15) < 1e-6) { ++loopValue; break;}
			++loopValue;
		}

		double currentBinValue = 0.15; //0.15 corresponds to 0.2 MeV bin
		 int currentBinNumber = 2;

		//Re-bin from 0.15 MeV - 0.95 MeV (0.2 - 0.9 MeV bins)
		while (std::getline(spaceSplitter,word, ' '))
		{
			double valueOver10 = std::stod(word)/10;

			Uwes_KE_spectrum_averaged.Fill(currentBinValue,valueOver10);

			// std::cout << "Current bin: " << currentBinValue+0.05 << " Value: " << word << " Averaged: " << valueOver10 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(currentBinNumber) << std::endl;

			if (innerLoopValue == 9) {innerLoopValue = -1; currentBinValue += 0.1; ++currentBinNumber; }
			if (std::fabs(uwesBinValues[loopValue] - 0.95) < 1e-6) { ++loopValue; innerLoopValue = 0; break;}
			++loopValue; ++innerLoopValue;
		}

		currentBinValue = 0.95; //Corresponds to 1 MeV bin

		//Re-bin from 0.95 MeV - 1.5 MeV (1 MeV bin)
		while (std::getline(spaceSplitter,word, ' '))
		{
			double valueOver55 = std::stod(word)/55;

			Uwes_KE_spectrum_averaged.Fill(currentBinValue,valueOver55);

			// std::cout << "Current bin: " << currentBinValue+0.05 << " Value: " << word << " Averaged: " << valueOver55 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(currentBinNumber) << std::endl;

			if (std::fabs(uwesBinValues[loopValue] - 1.5) < 1e-6) { ++loopValue; innerLoopValue = 0; ++currentBinNumber; break;}
			++loopValue; ++innerLoopValue;
		}

		currentBinValue = 1.5;

		//Re-bin from 1.5 MeV - 80 MeV (2 - 80 MeV bins)
		while (std::getline(spaceSplitter,word, ' '))
		{
			double valueOver100 = std::stod(word)/100;

			Uwes_KE_spectrum_averaged.Fill(currentBinValue,valueOver100);

			// std::cout << "Current bin: " << currentBinValue+0.5 << " Value: " << word << " Averaged: " << valueOver100 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(currentBinNumber) << std::endl;

			if (innerLoopValue == 99) {innerLoopValue = -1; currentBinValue += 1; ++currentBinNumber; }
			if (std::fabs(uwesBinValues[loopValue] - 80.5) < 1e-6) { ++loopValue; innerLoopValue = 0; break;}
			++loopValue; ++innerLoopValue;
		}

		//Keep running up until end of line 
		while (std::getline(spaceSplitter,word, ' '))
		{
			if (std::fabs(uwesBinValues[loopValue] - 100) < 1e-6) { ++loopValue; innerLoopValue = 0; break;}
			++loopValue;
		}

		spectra.push_back(Uwes_KE_spectrum_averaged);
	}

	return std::make_pair(distance,spectra);
}

void PlotUwesProtonSpectra(std::pair<std::vector<double>,std::vector<TH1F>>& UwesSpectra)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	double distance = 0.05;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (auto& protonSpectra:std::get<1>(UwesSpectra))
	{
		//Parse the output name
	
		std::stringstream nameStream;
		nameStream << std::setprecision(3) << distance;
		std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/KE_spectra/Uwe/" + nameStream.str() + "mm.jpg";

		//Set axes
		protonSpectra.SetTitle("");
		protonSpectra.GetYaxis()->SetTitle("frequency");
		protonSpectra.GetXaxis()->SetTitle("E [MeV]");
		protonSpectra.GetXaxis()->CenterTitle(true);
		protonSpectra.GetYaxis()->CenterTitle(true);
		protonSpectra.GetXaxis()->SetTitleFont(42);
		protonSpectra.GetYaxis()->SetTitleFont(42);
		protonSpectra.GetXaxis()->SetTitleSize(0.042);
		protonSpectra.GetYaxis()->SetTitleSize(0.048);
		protonSpectra.GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap
		// gPad->SetLogx();

		//Set histogram fill
		protonSpectra.SetFillColorAlpha(kAzure+3, 0.5);
		protonSpectra.SetLineColor(kBlack);
		protonSpectra.SetLineWidth(2);
		protonSpectra.SetLineStyle(1);

		//Draw the histogram
		protonSpectra.Draw("HIST");

		//Draw the inline title
		TLatex *t = new TLatex(0.015, 0.935, (TString)(nameStream.str() + " mm"));
		t->SetNDC(); //set position in coordinate system relative to canvas
		t->Draw();

		//Save
		c->SaveAs((TString)outputName);

		delete t; //Don't want to leak memory

		distance += 0.1;
	}
	delete c;
}

void Uwes_Spectrum()
{
	auto UwesSpectrumLibrary = Import_Spectrum_and_Rebin();
	PlotUwesProtonSpectra(UwesSpectrumLibrary);
}
