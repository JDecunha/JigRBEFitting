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

template<typename T>
class Logspace {
private:
    T curValue, base;

public:
    Logspace(T first, T base) : curValue(first), base(base) {}

    T operator()() {
        T retval = curValue;
        curValue *= base;
        return retval;
    }
};

std::vector<double> pyLogspace(double start, double stop, int num = 50, double base = 10) 
{
    double realStart = pow(base, start);
    double realBase = pow(base, (stop-start)/num);

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart,realBase));
    return retval;
}

TH1D Load_Martas_Measured_Spectrum(std::string filePath)
{
	//The path, an input stream, and a temp to hold the line being readout
	std::ifstream spectrumFile(filePath);
	std::string ReadoutLine;

	//Load bins
	//Loop counter
	int i = 0; 
	std::getline(spectrumFile,ReadoutLine); //Skip the header line
	std::vector<double> bins;

	while (std::getline(spectrumFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		int loopValue = 0;
		float energy = 0.;

		while (std::getline(spaceSplitter,word, ','))
		{
			if (word != "") //if it's not an empty space
			{
				++loopValue; //if it's not an empty space indicate we have gone on to the next word
				if (loopValue == 1) {bins.push_back(std::stod(word));break;}  //grab the bin
			}
		}
		
		//if(i > 10) {break;}
		++i;
	}

	//for (double en:bins) {std::cout << "bins: " << en << std::endl;}

	TH1D output = TH1D("Martas measured","Martas measured", 98, bins.data());

	//Load the data!
	spectrumFile.clear();
	spectrumFile.seekg(0);

	//Loop counter
	i = 1; 
	std::getline(spectrumFile,ReadoutLine); //Skip the header line

	while (std::getline(spectrumFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		int loopValue = 0;
		float energy = 0.;

		while (std::getline(spaceSplitter,word, ','))
		{
			if (word != "") //if it's not an empty space
			{
				++loopValue; //if it's not an empty space indicate we have gone on to the next word
				if (loopValue == 5) {output.SetBinContent(i,std::stod(word));break;}  //grab the bin
			}
		}
		
		//if(i > 10) {break;}
		++i;
	}

	//Utilities::PrintHistogram(output);
	//Normalizing the spectrum (so it adds to 1)
	double normalization = 0;

	for (int i = 1; i <= output.GetNbinsX(); ++i) //calculate the normalization value
	{
		normalization += output.GetBinContent(i)*output.GetBinWidth(i);
	}

	for (int i = 1; i <= output.GetNbinsX(); ++i) //normalize
	{
		double normalizedVal = output.GetBinContent(i)/normalization;
		output.SetBinContent(i,normalizedVal);
		//std::cout << KE_spectrum_averaged.GetBinCenter(i) << " " << KE_spectrum_averaged.GetBinContent(i) << std::endl;
	}

	return output;
}

TH1D Import_Spectrum(std::string filePath)
{
	//Let's just explicitly grab the energies from the library for Erik
	std::vector<double> binValues;

	//Hardcoded path and target size
	std::string fyFolder = "/home/joseph/Dropbox/Documents/Work/PHD_local/December2023_5um_logarithmic";
	std::string targetSize = "1e3"; //1 um diameter target spheres	

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	fyFolder = fyFolder + "/" + targetSize + "nm"; 

	for (const auto &entry : std::filesystem::directory_iterator(fyFolder)) //Loop over all the files in the folder
	{
		std::string energy = Utilities::GetFileEnergy(entry.path().filename()); //Get the energy of the file
		binValues.push_back(std::stod(energy));
	}
	std::sort(binValues.begin(), binValues.end());
	//Just to check the logspace worked correctly
	//for (float en:binValues) {std::cout << "Bin energy: " << en << std::endl;}

	//
	// Create the TH1F
	//

	TH1D Marta_KE_spectrum = TH1D("Ke spectrum","KE spectrum from Marta", 299, binValues.data());

	//The path, an input stream, and a temp to hold the line being readout
	std::ifstream spectrumFile(filePath);
	std::string ReadoutLine;

	//Loop counter
	int i = 0; 

	while (std::getline(spectrumFile,ReadoutLine)) //Split the file by lines
	{
		//Split the line by words
		std::istringstream spaceSplitter(ReadoutLine); std::string word;
		int loopValue = 0;
		float energy = 0.;

		while (std::getline(spaceSplitter,word, ' '))
		{
			if (word != "") //if it's not an empty space
			{
				++loopValue; //if it's not an empty space indicate we have gone on to the next word
				if (loopValue == 6) {energy = std::stod(word);}  //grab the energy
				if (loopValue == 8) 
				{
					if (word == "2212") //if it's a proton add to the list
					{
						Marta_KE_spectrum.Fill(energy);
					} 

				break;} 
			}
		}
		
		//if(i > 10) {break;}
		++i;
	}

	//Just to check the energies were loaded correctly
	//for (float en:protonEnergies) {std::cout << "energy: " << en << std::endl;}

	//Marta_KE_spectrum.SetBinContent(10,42);
	//Utilities::PrintHistogram(Marta_KE_spectrum);

	///
	/// Last thing to do is just to normalize the spectrum
	///

	//Normalizing the spectrum (so it adds to 1)
	double normalization = 0;

	//Do we want it to be a PDF or a count spectrum?
	/*for (int i = 1; i <= Marta_KE_spectrum.GetNbinsX(); ++i) //calculate the normalization value
	{
		auto val = Marta_KE_spectrum.GetBinContent(i);
		auto width = Marta_KE_spectrum.GetBinWidth(i);
		Marta_KE_spectrum.SetBinContent(i,(val/width));
	}*/

	for (int i = 1; i <= Marta_KE_spectrum.GetNbinsX(); ++i) //calculate the normalization value
	{
		normalization += Marta_KE_spectrum.GetBinContent(i);
	}

	for (int i = 1; i <= Marta_KE_spectrum.GetNbinsX(); ++i) //normalize
	{
		double normalizedVal = Marta_KE_spectrum.GetBinContent(i)/normalization;
		Marta_KE_spectrum.SetBinContent(i,normalizedVal);
		//std::cout << KE_spectrum_averaged.GetBinCenter(i) << " " << KE_spectrum_averaged.GetBinContent(i) << std::endl;
	}

	//Utilities::PrintHistogram(Marta_KE_spectrum);

	return Marta_KE_spectrum;
}

void PlotProtonSpectra(TH1D protonSpectra, std::string name)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	// double distance = 0.05;

	//Parse the output name
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/KE_spectra/Marta/" + name + ".jpg";

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
	//TLatex *t = new TLatex(0.015, 0.935, (TString)(nameStream.str() + " mm"));
	//t->SetNDC(); //set position in coordinate system relative to canvas
	//t->Draw();

	//Save
	c->SaveAs((TString)outputName);

	//delete t; //Don't want to leak memory


	delete c;
}

TH1D GetKeWeightedFrequencyLinealSpectra(TH1D KESpectrum)
{

	//Specify hardcoded path to our f(y) library and target size
	// std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string fyFolder =  "/home/joseph/Dropbox/Documents/Work/PHD_local/December2023_5um_logarithmic";

	std::string targetSize = "1e3";

	//Just get the bins for the output spectrum
	TH1D outputLinealSpectrum = utils::GetFy(fyFolder, 0.001, targetSize);
	outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

	//Loop over the KE spectrum and weight the f(y) spectra accordingly
	for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 1 to skip underflow
	{
		//Because of the way we re-defined the KE histogram, high edge actually corresponds to the middle
		double lowEdge = KESpectrum.GetBinLowEdge(i);// + KESpectrum.GetBinWidth(i);
		//std::cout << lowEdge << std::endl;
		TH1D tempLinealSpectrum = utils::GetNy(fyFolder, lowEdge, targetSize); //Get the N(y) spectrum at that energy

		//Scale the spectrum by the fluence
		tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

		//Scale the spectrum so N(y) matches an effective count of 10 million
		long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, lowEdge, targetSize);
		double scalingFactor = double(1e7)/double(effectiveNumTracks);
		tempLinealSpectrum.Scale(scalingFactor);

		//Add the weighted spectrum to the output fy spectrum
		outputLinealSpectrum.Add(&tempLinealSpectrum);
	}	
	
	//Convert N(y) to f(y)
	utils::PMF_to_FrequencyFunction(&outputLinealSpectrum);
	//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7

	return outputLinealSpectrum;
}

TH1D GetKeWeightedFrequencyLinealSpectraEnergyThreshold(TH1D KESpectrum, double EnergyThreshold)
{

	//Specify hardcoded path to our f(y) library and target size
	// std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string fyFolder =  "/home/joseph/Dropbox/Documents/Work/PHD_local/December2023_5um_logarithmic";

	std::string targetSize = "1e3";

	//Just get the bins for the output spectrum
	TH1D outputLinealSpectrum = utils::GetFy(fyFolder, 0.001, targetSize);
	outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

	//Loop over the KE spectrum and weight the f(y) spectra accordingly
	for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 1 to skip underflow
	{
		//Because of the way we re-defined the KE histogram, high edge actually corresponds to the middle
		double lowEdge = KESpectrum.GetBinLowEdge(i);// + KESpectrum.GetBinWidth(i);
		//std::cout << lowEdge << std::endl;
		TH1D tempLinealSpectrum = utils::GetNy(fyFolder, lowEdge, targetSize); //Get the N(y) spectrum at that energy

		//Scale the spectrum by the fluence
		tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

		//Scale the spectrum so N(y) matches an effective count of 10 million
		long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, lowEdge, targetSize);
		double scalingFactor = double(1e7)/double(effectiveNumTracks);
		tempLinealSpectrum.Scale(scalingFactor);

		//Add the weighted spectrum to the output fy spectrum
		outputLinealSpectrum.Add(&tempLinealSpectrum);
	}	

	// Energy threshold
	int length = outputLinealSpectrum.GetNbinsX();
	for (int i = 1; i <= length; i++) //bins per bin width
	{
		double low_edge = outputLinealSpectrum.GetBinLowEdge(i);

		if (low_edge < EnergyThreshold)
		{
			outputLinealSpectrum.SetBinContent(i,0);
		}

	}
	
	//Convert N(y) to f(y)
	utils::PMF_to_FrequencyFunction(&outputLinealSpectrum);
	//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7

	return outputLinealSpectrum;
}


void  PlotKeWeightedLinealSpectra(TH1D spectra, std::string name)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right

	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);


	//Get the d(y) spectrum and multiply by y
	//TH1D* spectra = new TH1D(spectrum);
	// utils::PMF_to_DoseFunction(spectra);
	//utils::PrintHistogram(spectra);
	utils::Prepare_for_Semilog(&spectra); //Multiply by y, to get y*f(y)

	// std::string name = std::to_string(distance);
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/" + name + ".jpg";

	//Draw 
	spectra.Draw("HIST");

	//Set axis settings
	spectra.SetTitle("");
	spectra.SetTitleSize(0.03,"t"); //this doesn't do anything
	spectra.GetYaxis()->SetTitle("y #upoint f(y)");
	spectra.GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
	spectra.GetXaxis()->CenterTitle(true);
	spectra.GetYaxis()->CenterTitle(true);
	spectra.GetXaxis()->SetTitleFont(42);
	spectra.GetYaxis()->SetTitleFont(42);
	spectra.GetXaxis()->SetTitleSize(0.042);
	spectra.GetYaxis()->SetTitleSize(0.048);
	spectra.GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap
	gPad->SetLogx();

	//Set the spectrum fill and line settings
	spectra.SetFillColorAlpha(kAzure+3, 0.5);
	spectra.SetLineColor(kBlack);
	spectra.SetLineWidth(2);
	spectra.SetLineStyle(1);

	//Save
	c->SaveAs((TString)outputName);

	//Delete so we don't leak memory
	//delete spectra;
	delete c;
}

void Martas_vs_Calculated(TH1D* h, TH1D* h2, std::string outputLocation)
{
		gStyle->SetOptStat(0); //Don't print the stats window in the top right

		TCanvas* c = new TCanvas("c","c");
		c->SetCanvasSize(2040, 1640);
		c->SetWindowSize(2040./2, 1640./2);
		c->SetLeftMargin(0.15);
		c->SetBottomMargin(0.15);

		utils::VerifyNormalization(*h); 
		utils::VerifyNormalization(*h2); //Shows all spectra are normalized to 1 within 1e-7

		Utilities::Prepare_for_Semilog(h);
		Utilities::Prepare_for_Semilog(h2);

		gStyle->SetTitleFont(42,"t");

		h->Draw("hist");
		h2->Draw("same hist");
		
		//Set range
		//h2->GetYaxis()->SetRangeUser(0,1);
		// h2->GetXaxis()->SetRangeUser(10,150);
		h->GetYaxis()->SetRangeUser(0,0.6);

		//Set titles
		//std::string title = std::to_string(int(energy));
		//title += " MeV";
		h->SetTitle("");
		h->SetFillColorAlpha(kAzure+3, 0.5);

		h->SetTitle("");
		h->SetTitleSize(0.03,"t"); //this doesn't do anything
		h->GetYaxis()->SetTitle("y #upoint f(y)");
		h->GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
		h->GetXaxis()->CenterTitle(true);
		h->GetYaxis()->CenterTitle(true);
		h->GetXaxis()->SetTitleFont(42);
		h->GetYaxis()->SetTitleFont(42);
		h->GetXaxis()->SetTitleSize(0.042);
		h->GetYaxis()->SetTitleSize(0.048);
		h->GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap
		//h2->GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap

		//Center
		h->GetXaxis()->CenterTitle(true);
		h->GetYaxis()->CenterTitle(true);

		//Offset x axis so no overlap
		h->GetXaxis()->SetTitleOffset(1.1);
		h->GetYaxis()->SetTitleOffset(1.0);


		//h->SetFillColorAlpha(kRed, 0.4);
		//h->SetFillStyle(3001);
		h->SetLineColor(kBlue);
		h->SetLineWidth(0);

		//h2->SetFillColorAlpha(kAzure+3, 0.5);
		//h->SetFillStyle(3001);
		h2->SetLineColor(kBlack);
		h2->SetLineWidth(5);
		// h->SetLineWidth(5);
		// h->SetLineStyle(9);

		gPad->SetLogx();

		double xstart = 0.64;
		double xwidth = 0.24;
		double ystart = 0.725;
		double ywidth = 0.15;
		auto legend = new TLegend(xstart,ystart,xstart+xwidth,ystart+ywidth);//x start, y start, x end, yend
		// auto legend = new TLegend(0.70,0.85,0.70+0.25,0.85+0.13);
		//legend->SetHeader("","C"); // option "C" allows to center the header
		legend->AddEntry(h,"Summation","f");
		legend->AddEntry(h2,"Measurement","L");
		legend->SetTextSize(0.035);
		legend->Draw();

		c->SaveAs((TString)outputLocation);
}

void Martas_InField()
{
	// 3.3 cm

	std::string FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D33_L00/PhaseSpace.phsp";
	std::string Outname = "33mmInField";
	std::string MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/33mm/f01.csv";
	
	auto protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	auto fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	auto Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");


	// 6.6 cm

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D6.0_L0.0/PhaseSpace.phsp";
	Outname = "66mmInField";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/66mm/histogram_f02.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 11.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D10.7_L0.0/PhaseSpace.phsp";
	Outname = "116mmInField";
	MartasData = "//home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/116mm/histogram_f07.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 12.9(?) or 12.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D12.0_L0.0/PhaseSpace.phsp";
	Outname = "129mmInField";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/129or126mm/histogram_f11.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

}

void Martas_InField_EnergyThreshold()
{
	// 3.3 cm

	std::string FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D33_L00/PhaseSpace.phsp";
	std::string Outname = "33mmInField_Threshold";
	std::string MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/33mm/f01.csv";
	
	auto protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	auto fy = GetKeWeightedFrequencyLinealSpectraEnergyThreshold(protonSpectrum,0.1);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	auto Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");


	// 6.6 cm

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D6.0_L0.0/PhaseSpace.phsp";
	Outname = "66mmInField_Threshold";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/66mm/histogram_f02.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectraEnergyThreshold(protonSpectrum,0.1);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 11.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D10.7_L0.0/PhaseSpace.phsp";
	Outname = "116mmInField_Threshold";
	MartasData = "//home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/116mm/histogram_f07.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectraEnergyThreshold(protonSpectrum,0.1);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 12.9(?) or 12.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D12.0_L0.0/PhaseSpace.phsp";
	Outname = "129mmInField_Threshold";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/129or126mm/histogram_f11.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectraEnergyThreshold(protonSpectrum,0.1);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

}


void Martas_CubicTEPC()
{
	// 3.3 cm

	std::string FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/CubicTEPC/D2.7_L0/PhaseSpace.phsp";
	std::string Outname = "33mmInField_CubicTEPC";
	std::string MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/33mm/f01.csv";
	
	auto protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	auto fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	auto Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 11.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/CubicTEPC/D10.7_L0/PhaseSpace.phsp";
	Outname = "116mmInField_CubicTEPC";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/116mm/histogram_f07.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

	// 12.9(?) or 12.6 cm 

	FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/CubicTEPC/D12.0_L0/PhaseSpace.phsp";
	Outname = "129mmInField_CubicTEPC";
	MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/129or126mm/histogram_f11.csv";

	protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");

}

void Martas_OutOfField()
{
	// 11.6 cm out of field

	std::string FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D11.0_L5.0/PhaseSpace.phsp";
	std::string Outname = "116mm_5mmlateral_OutOfField";
	std::string MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/116mm5mmlateral/histogram_f16.csv";
	
	auto protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	auto fy = GetKeWeightedFrequencyLinealSpectra(protonSpectrum);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	auto Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");
}

void Martas_OutOfField_EnergyThreshold()
{
	// 11.6 cm out of field

	std::string FilePath = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/Phasespaces/RegularTEPC/D11.0_L5.0/PhaseSpace.phsp";
	std::string Outname = "116mm_5mmlateral_OutOfField_EnergyThreshold";
	std::string MartasData = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/Data from Marta/Combined/TEPC_Measurement/116mm5mmlateral/histogram_f16.csv";
	
	auto protonSpectrum = Import_Spectrum(FilePath);
	PlotProtonSpectra(protonSpectrum, Outname);

	auto fy = GetKeWeightedFrequencyLinealSpectraEnergyThreshold(protonSpectrum, 0.1);
	PlotKeWeightedLinealSpectra(fy, Outname);
	
	auto Measured = Load_Martas_Measured_Spectrum(MartasData);
	
	Martas_vs_Calculated(&fy,&Measured,"/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/Marta/"+Outname+"_comparison.jpg");
}




void Martas_Spectrum()
{
	Martas_InField_EnergyThreshold();
}
