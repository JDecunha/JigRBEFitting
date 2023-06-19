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
#include "TGraph.h"
//This project
#include "Uwes_Spectrum.h"
#include "Utilities.h"
#include "BWF_Fitting_Results.h"
#include "LETCalculator.h"

namespace utils = Utilities;

//Alright, Uwe's spectrum
//10 keV energy resolution. From 0.01 MeV up to 100 MeV (10 k bins)
//1 mm spatial resolution. From 1 mm to 70 mm (70 spectra)

std::vector<std::pair<float,TH1F>> Import_Spectrum_and_Rebin()
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
	std::vector<double> distance = utils::Linspace(0.5,69.5, 70);
	std::reverse(distance.begin(), distance.end());
	TH1::AddDirectory(false); 
	std::vector<std::pair<float,TH1F>> spectra;

	//Loop counter
	int i = 0; 

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
			// if (loopValue == 0) {std::cout << "first bin val: " << word << std::endl;}
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

		//Normalizing the spectrum (so it adds to 1)
		float normalization = 0;

		for (int i = 0; i <= Uwes_KE_spectrum_averaged.GetNbinsX(); ++i) //calculate the normalization value
		{
			normalization += Uwes_KE_spectrum_averaged.GetBinContent(i);
		}

		for (int i = 0; i <= Uwes_KE_spectrum_averaged.GetNbinsX(); ++i) //normalize
		{
			float normalizedVal = Uwes_KE_spectrum_averaged.GetBinContent(i)/normalization;
			Uwes_KE_spectrum_averaged.SetBinContent(i,normalizedVal);
			//std::cout << KE_spectrum_averaged.GetBinCenter(i) << " " << KE_spectrum_averaged.GetBinContent(i) << std::endl;
		}

			auto outpair = std::make_pair(distance[i], Uwes_KE_spectrum_averaged);
			spectra.push_back(outpair);

			++i;
		}

	return spectra;
}

std::vector<std::pair<float,TH1F>> Import_Spectrum_and_Rebin_01MeV()
{
	//The path, an input stream, and a temp to hold the line being readout
	std::string filePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Documents/Uwes_First_Spectrum.dat";
	std::ifstream spectrumFile(filePath);
	std::string ReadoutLine;

	//Create doubles corresponding to Uwe's bins
	std::vector<double> uwesBinValues = utils::Linspace(0.01,100, 10000);

	//Create the newBins we want to output with
	std::vector<double> newBins = utils::Linspace(0,80,801);

	//The output of this function
	std::vector<double> distance = utils::Linspace(0.5,69.5, 70);
	std::reverse(distance.begin(), distance.end());
	TH1::AddDirectory(false); 
	std::vector<std::pair<float,TH1F>> spectra;

	//Loop counter
	int i = 0; 

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
			// if (loopValue == 0) {std::cout << "first bin val: " << word << std::endl;}
			double valueOver15 = std::stod(word)/15;
			Uwes_KE_spectrum_averaged.Fill(0.05,valueOver15); //0.05 corresponds to 0.1 MeV bin
			// std::cout << "Value: " << word << " Averaged: " << valueOver15 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(1) << std::endl;
			// std::cout << uwesBinValues[loopValue] << std::endl;
			if (std::fabs(uwesBinValues[loopValue] - 0.15) < 1e-6) { ++loopValue; break;}
			++loopValue;
		}

		double currentBinValue = 0.15; //0.15 corresponds to 0.2 MeV bin
		int currentBinNumber = 2;

		//Re-bin from 0.15 MeV - 80 MeV MeV 
		while (std::getline(spaceSplitter,word, ' '))
		{
			double valueOver10 = std::stod(word)/10;

			Uwes_KE_spectrum_averaged.Fill(currentBinValue,valueOver10);

			// std::cout << "Current bin: " << currentBinValue+0.05 << " Value: " << word << " Averaged: " << valueOver10 << " Current bin value: " << Uwes_KE_spectrum_averaged.GetBinContent(currentBinNumber) << std::endl;

			if (innerLoopValue == 9) {innerLoopValue = -1; currentBinValue += 0.1; ++currentBinNumber; }
			if (std::fabs(uwesBinValues[loopValue] - 80.5) < 1e-6) { ++loopValue; innerLoopValue = 0; break;}
			++loopValue; ++innerLoopValue;
		}

		//Keep running up until end of line 
		while (std::getline(spaceSplitter,word, ' '))
		{
			if (std::fabs(uwesBinValues[loopValue] - 100) < 1e-6) { ++loopValue; innerLoopValue = 0; break;}
			++loopValue;
		}

		//Normalizing the spectrum (so it adds to 1)
		float normalization = 0;

		for (int i = 0; i <= Uwes_KE_spectrum_averaged.GetNbinsX(); ++i) //calculate the normalization value
		{
			normalization += Uwes_KE_spectrum_averaged.GetBinContent(i);
		}

		for (int i = 0; i <= Uwes_KE_spectrum_averaged.GetNbinsX(); ++i) //normalize
		{
			float normalizedVal = Uwes_KE_spectrum_averaged.GetBinContent(i)/normalization;
			Uwes_KE_spectrum_averaged.SetBinContent(i,normalizedVal);
			//std::cout << KE_spectrum_averaged.GetBinCenter(i) << " " << KE_spectrum_averaged.GetBinContent(i) << std::endl;
		}

			auto outpair = std::make_pair(distance[i], Uwes_KE_spectrum_averaged);
			spectra.push_back(outpair);

			++i;
		}

	return spectra;
}

void  SaveUweSpectra()
{
	gInterpreter->GenerateDictionary("pair<float,TH1F>;vector<pair<float,TH1F> >","TH1.h;float;utility;vector");
	auto spectra = Import_Spectrum_and_Rebin();
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/UwesFirstProtonSpectra.root","RECREATE");
	spectrumFile->WriteObject(&spectra, "Spectrum_Library");
	delete spectrumFile;
}

void PlotUwesProtonSpectra(std::vector<std::pair<float,TH1F>>& UwesSpectra)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	// double distance = 0.05;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (auto& spectrumAtDistance:UwesSpectra)
	{
		//Parse the output name
		double distance = std::get<0>(spectrumAtDistance);
		TH1F protonSpectra = std::get<1>(spectrumAtDistance);
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

std::vector<std::pair<float,TH1D>> GetUwesKeWeightedFrequencyLinealSpectra(std::string targetSize)
{
	/*
	Retrieves the lineal energy spectra from Uwes first spectrum
	*/

	//Specify hardcoded path to our f(y) library and target size
	// std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string fyFolder = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/5um_0to100MeV_April2023Library";

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<float,TH1F>;vector<pair<float,TH1F> >","TH1.h;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<float,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/UwesFirstProtonSpectra.root");
	spectrumFile->GetObject("Spectrum_Library",keSpectra);
	spectrumFile->Close();
	
	//This will hold the output
	std::vector<std::pair<float,TH1D>> linealEnergyLibrary; 

	
	for(const auto& spectra:*keSpectra)
	{
		//Step 1: Weight the lineal energy spectrum by the KE spectrum.
		//Grab the KE specta for a given irradiation
		TH1F KESpectrum = std::get<1>(spectra);
		double KEDistance = std::get<0>(spectra);

		//Just get the bins for the output spectrum
		TH1D outputLinealSpectrum = utils::GetFy(fyFolder, 0.1, targetSize);
		outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

		//Loop over the KE spectrum and weight the f(y) spectra accordingly
		for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 1 to skip underflow
		{
			//Because of the way we re-defined the KE histogram, high edge actually corresponds to the middle
			double highEdge = KESpectrum.GetBinLowEdge(i) + KESpectrum.GetBinWidth(i);
			TH1D tempLinealSpectrum = utils::GetNy(fyFolder, highEdge, targetSize); //Get the N(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Scale the spectrum so N(y) matches an effective count of 10 million
			long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, highEdge, targetSize);
			double scalingFactor = double(1e7)/double(effectiveNumTracks);
			tempLinealSpectrum.Scale(scalingFactor);

			//Add the weighted spectrum to the output fy spectrum
			outputLinealSpectrum.Add(&tempLinealSpectrum);
		}	
		
		//Convert N(y) to f(y)
		utils::PMF_to_FrequencyFunction(&outputLinealSpectrum);
		//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7
		linealEnergyLibrary.push_back(std::make_pair<float,TH1D>(std::move(KEDistance),std::move(outputLinealSpectrum)));
	}

	delete keSpectra; //we own the keSpectra pointer so we have to delete it

	return linealEnergyLibrary;
}

void  PlotKeWeightedLinealSpectra(const std::vector<std::pair<float,TH1D>>& linealLibrary)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right

	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (const auto& linealSpectra:linealLibrary)
	{
			//Get the d(y) spectrum and multiply by y
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			// utils::PMF_to_DoseFunction(spectra);
			utils::Prepare_for_Semilog(spectra); //Multiply by y, to get y*f(y)

			//Get the name corresponding to the jig position
			double distance = std::get<0>(linealSpectra);
			std::string s;
		    std::stringstream sstream;
		    sstream.setf(std::ios::fixed);
		    sstream.precision(1);
		    sstream << distance;
		    std::string name = sstream.str();
			// std::string name = std::to_string(distance);
			std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Uwes/" + name + "mm.jpg";

			//Draw 
			spectra->Draw("HIST");

			//Set axis settings
			spectra->SetTitle("");
			spectra->SetTitleSize(0.03,"t"); //this doesn't do anything
			spectra->GetYaxis()->SetTitle("y #upoint f(y)");
			spectra->GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
			spectra->GetXaxis()->CenterTitle(true);
			spectra->GetYaxis()->CenterTitle(true);
			spectra->GetXaxis()->SetTitleFont(42);
			spectra->GetYaxis()->SetTitleFont(42);
			spectra->GetXaxis()->SetTitleSize(0.042);
			spectra->GetYaxis()->SetTitleSize(0.048);
			spectra->GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap
			gPad->SetLogx();

			//Set the spectrum fill and line settings
			spectra->SetFillColorAlpha(kAzure+3, 0.5);
			spectra->SetLineColor(kBlack);
			spectra->SetLineWidth(2);
			spectra->SetLineStyle(1);

			//Draw the inline title
			name = name + " mm";
			TLatex *t = new TLatex(0.015, 0.935, (TString)name);
			t->SetNDC(); //set position in coordinate system relative to canvas
			t->Draw();

			//Save
			c->SaveAs((TString)outputName);

			//Delete so we don't leak memory
			delete spectra;
			delete t;
	}

	delete c;
}

void Plot_RBE_SF10_Uwes_Spectrum(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, std::vector<std::pair<float,TH1D>>& linealLibrary, double const* alphaBetasCesium, BWF_Fitting_Results results)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	int i = 0;

	std::vector<float> xAxis;
	std::vector<float> toPlot;
	std::vector<double> alphaBetasProton;

	//Iterate through each of the experiments levels
	for(const std::pair<float,TH1D>& spectrumPair:linealLibrary) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaValue = 0; double betaValue = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryAlphaVal = results.alphaFunc.GetValue(center);
			alphaValue += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = results.betaFunc.GetValue(center);
			betaValue += ryBetaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)
		}
		
		alphaBetasProton.push_back(alphaValue);
		alphaBetasProton.push_back(betaValue);

		++i;
	}

	for (int i = 0; i < linealLibrary.size(); ++i)
	{
		//quadratic formula
		double b = alphaBetasCesium[0];
		double a = alphaBetasCesium[1];
		double c = -1;
		double resultCesium = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		b = alphaBetasProton[i*2];
		a = alphaBetasProton[(i*2)+1];
		double resultProton = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Push back the ratio of doses
		toPlot.push_back(resultCesium/resultProton);
		xAxis.push_back(std::move(std::get<0>(linealLibrary[i])));
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(linealLibrary.size(),xAxis.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("d [mm]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0.65);
	// gPad->SetLogy(); //set the new y Limits
	

	gr->GetXaxis()->SetRangeUser(0,48); //set the new x Limits
	gr->GetYaxis()->SetRangeUser(0,5.5); //set the new y Limits
	

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void Plot_RBE_SF10_Uwes_Spectrum_LET(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, std::vector<std::pair<float,float>>& letLibrary, double const* alphaBetasCesium, BWF_Fitting_Results results)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	std::vector<float> xAxis;
	std::vector<float> toPlot;
	std::vector<double> alphaBetasProton;

	//Iterate through each of the experiments levels
	for(const std::pair<float,float>& letPair:letLibrary) //First we iterate over each of the lineal energy spectra
	{

		double LETd = std::get<1>(letPair); 
		double alphaValue =  results.alphaFunc.GetValue(LETd);
		double betaValue = results.betaFunc.GetValue(LETd);
		
		alphaBetasProton.push_back(alphaValue);
		alphaBetasProton.push_back(betaValue);
	}

	for (int i = 0; i < letLibrary.size(); ++i)
	{
		//quadratic formula
		double b = alphaBetasCesium[0];
		double a = alphaBetasCesium[1];
		double c = -1;
		double resultCesium = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		b = alphaBetasProton[i*2];
		a = alphaBetasProton[(i*2)+1];
		double resultProton = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Push back the ratio of doses
		toPlot.push_back(resultCesium/resultProton);
		xAxis.push_back(std::move(std::get<0>(letLibrary[i])));
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(letLibrary.size(),xAxis.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("d [mm]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(-1);
	gr->GetXaxis()->SetTitleOffset(1.15);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void Uwes_Spectrum()
{
	// auto UwesSpectrumLibrary = Import_Spectrum_and_Rebin_01MeV();
	// PlotUwesProtonSpectra(UwesSpectrumLibrary);
	auto uwesLib = GetUwesKeWeightedFrequencyLinealSpectra("1e3");
	// PlotKeWeightedLinealSpectra(uwesLib);

	double* cesiumAlphaBeta = new double[2];
	cesiumAlphaBeta[0] = 0.05; cesiumAlphaBeta[1] = 0.041;

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->SetBottomMargin(0.13);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	auto legend = new TLegend(0.14,0.76,0.14+0.23,0.62+0.23);
	legend->SetTextSize(0.04);
	TAttMarker markerAtts;
	TAttLine lineStyle{};
	markerAtts.SetMarkerSize(10);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	BWF_Fitting_Results CubicCubicFyBestResults;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);
	CubicBWF.SetValues(std::vector<double> {0.00013752, -0.00399419, 0.0252074, 0.0564586});
	CubicCubicFyBestResults.alphaFunc = CubicBWF;
	CubicBWF.SetValues(std::vector<double> {9.41183e-06, -0.000406718, 0.00310442, 0.0296854});
	CubicCubicFyBestResults.betaFunc = CubicBWF;

	markerAtts.SetMarkerColor(kViolet+4);
	markerAtts.SetMarkerStyle(20);
	Plot_RBE_SF10_Uwes_Spectrum(c, legend, "f(y) spectrum", markerAtts, "AP", uwesLib,cesiumAlphaBeta,CubicCubicFyBestResults);

	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	LE2BWF.SetValues(std::vector<double> {0.0941433, -0.00583604, 0.00233011});

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);
	FifthBWF.SetValues(std::vector<double> {1.10568e-06, -4.67063e-05, 0.000773506, -0.00602498, 0.0228525, 0.00195857});

	BWF_Fitting_Results BestLETH1437;
	BestLETH1437.alphaFunc = LE2BWF;
	BestLETH1437.betaFunc = FifthBWF;

	auto LETLibrary = ComputeUwesLETs();
	markerAtts.SetMarkerColor(kTeal+3);
	markerAtts.SetMarkerStyle(34);
	Plot_RBE_SF10_Uwes_Spectrum_LET(c, legend, "LET_{d}", markerAtts, "P", LETLibrary,cesiumAlphaBeta,BestLETH1437);

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/RBE_10_SF_Uwes_H1437_notitle.jpg";
	c->SaveAs((TString)outputName); 
}
