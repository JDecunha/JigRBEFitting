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

//This project
#include "LinealSpectra.h"
#include "Utilities.h"

namespace utils = Utilities; //namespace alias

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetMonoenergeticLinealSpectra()
{
	/*
	Retrieves the lineal energy spectra for 0.1 MeV - 79 MeV monoenergetic protons.
	These energies correspond to the bins in Fada's KE spectrum
	*/

	//Hardcoded path and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string targetSize = "1e3"; //1 um diameter target spheres	

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	std::vector<std::pair<std::string, TH1D>> monoenergeticSpectra; //List to hold our f(y) spectrum
	std::vector<std::pair<std::string,TH1F>>* keSpectra; //A pointer to Fada's KE spectrum

	//Pull the data from our KE-spectrum file.
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/FadasProtonSpectra.root");
	spectrumFile->GetObject("Spectrum_Library",keSpectra);
	spectrumFile->Close(); //Close file
	
	//Loop through all the bins of the KE spectrum
	for (int i = 2; i <= std::get<1>((*keSpectra)[0]).GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
	{
		//Because of the way we re-defined the KE histogram, "GetBinLowEdge" corresponds to the middle
		double energy = std::get<1>((*keSpectra)[0]).GetBinLowEdge(i);

		//utils::VerifyNormalization(utils::GetDy(fyFolder, energy, targetSize)); If you want to verify the normalization

		monoenergeticSpectra.push_back(std::make_pair<std::string,TH1D>(std::to_string(energy),utils::GetDy(fyFolder, energy, targetSize)));
	}

	delete keSpectra; //Since we own keSpectra we should delete it. std::vector's destructor will take care of std::strings and TH1Fs in it

	return monoenergeticSpectra;
}

void LinealSpectra::PlotMonoenergeticLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& monoenergeticSpectra)
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (const auto& linealSpectra:monoenergeticSpectra)
	{
		//Get the spectrum and multiply it by y, for semilog plotting
		TH1D spectra = std::get<1>(linealSpectra);
		utils::Prepare_for_Semilog(&spectra);

		//Parse the output name
		double energy = std::stod(std::get<0>(linealSpectra));
		std::stringstream nameStream;
		nameStream << std::setprecision(2) << energy;
		std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/monoSpectra/" + nameStream.str() + "MeV.jpg";

		//Set axes
		spectra.SetTitle("");
		spectra.GetYaxis()->SetTitle("y #upoint d(y)");
		spectra.GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
		spectra.GetXaxis()->CenterTitle(true);
		spectra.GetYaxis()->CenterTitle(true);
		spectra.GetXaxis()->SetTitleFont(42);
		spectra.GetYaxis()->SetTitleFont(42);
		spectra.GetXaxis()->SetTitleSize(0.042);
		spectra.GetYaxis()->SetTitleSize(0.048);
		spectra.GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap
		gPad->SetLogx();

		//Set histogram fill
		spectra.SetFillColorAlpha(kAzure+3, 0.5);
		spectra.SetLineColor(kBlack);
		spectra.SetLineWidth(2);
		spectra.SetLineStyle(1);

		//Draw the histogram
		spectra.Draw("HIST");

		//Draw the inline title
		TLatex *t = new TLatex(0.015, 0.935, (TString)(nameStream.str() + " MeV"));
		t->SetNDC(); //set position in coordinate system relative to canvas
		t->Draw();

		//Save
		c->SaveAs((TString)outputName);

		delete t; //Don't want to leak memory
	}
	delete c;
}

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetKeWeightedLinealSpectra()
{
	/*
	Retrieves the lineal energy spectra calculated at every location of the jig.
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string targetSize = "1e3"; //1 um diameter target spheres

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<std::string,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/FadasProtonSpectra.root");
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
		TH1D outputLinealSpectrum = utils::GetDy(fyFolder, 0.1, targetSize);
		outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

		//Loop over the KE spectrum and weight the f(y) spectra accordingly
		for (int i = 2; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			//Because of the way we re-defined the KE histogram, "GetBinLowEdge" actually corresponds to the middle
			TH1D tempLinealSpectrum = utils::GetDy(fyFolder, KESpectrum.GetBinLowEdge(i), targetSize); //Get the f(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Add the weighted spectrum to the output fy spectrum
			outputLinealSpectrum.Add(&tempLinealSpectrum);
		}	
		
		//For verifying normalization
		//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7

		linealEnergyLibrary.push_back(std::make_pair<std::string,TH1D>(std::move(KEColumnName),std::move(outputLinealSpectrum)));
	}

	delete keSpectra; //we own the keSpectra pointer so we have to delete it

	return linealEnergyLibrary;
}

void LinealSpectra::PlotKeWeightedLinealSpectra(const std::vector<std::pair<std::string,TH1D>>& linealLibrary)
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
			utils::Prepare_for_Semilog(spectra); //Multiply by y, to get y*d(y)

			//Get the name corresponding to the jig position
			std::string name = std::get<0>(linealSpectra);
			std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/" + name + ".jpg";

			//Draw 
			spectra->Draw("HIST");

			//Set axis settings
			spectra->SetTitle("");
			spectra->SetTitleSize(0.03,"t"); //this doesn't do anything
			spectra->GetYaxis()->SetTitle("y #upoint d(y)");
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

void LinealSpectra::PlotKeWeightedLinealSpectraMultigraph(const std::vector<std::pair<std::string,TH1D>>& linealLibrary)
{
	//FYI: this function leaks memory because the TH1Ds and TLatexes don't get deleted
	//But if you delete them in the loop they get deleted before the multigraph is printed
	//It's fine though because this is just plotting code.
	//Also this isn't a problem on the non-multigraph version, because we can delete them after plotting.

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->Divide(4,3);

	int i = 1;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (const auto& linealSpectra:linealLibrary)
	{
			c->cd(i); ++i;//Change the window of the canvas we're drawing on (since it's a multigraph)

			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			utils::Prepare_for_Semilog(spectra); //Multiply by y, to get y*d(y)
			std::string name = std::get<0>(linealSpectra);

			//Draw
			spectra->Draw("HIST");

			//Set axes
			spectra->SetTitle("");
			spectra->SetTitleSize(0.03,"t"); //this doesn't do anything
			spectra->GetYaxis()->SetTitle("y #upoint d(y)");
			spectra->GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
			spectra->GetXaxis()->CenterTitle(true);
			spectra->GetYaxis()->CenterTitle(true);
			spectra->GetXaxis()->SetTitleFont(42);
			spectra->GetYaxis()->SetTitleFont(42);
			spectra->GetXaxis()->SetTitleSize(0.042);
			spectra->GetYaxis()->SetTitleSize(0.048);
			gPad->SetLogx();

			//Set histogram fill and line settings
			spectra->SetFillColorAlpha(kAzure+3, 0.5);
			spectra->SetLineColor(kBlack);
			spectra->SetLineWidth(2);
			spectra->SetLineStyle(1);

			//Draw the inline title
			TLatex *t = new TLatex(0.015, 0.935, (TString)name);
			t->SetNDC(); //set position in coordinate system relative to canvas
			t->Draw();
	}

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/multigraphJigLinealSpectra.jpg";
	c->SaveAs((TString)outputName);

	delete c;
}

void LinealSpectra::SaveKeWeightedLinealSpectra()
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	auto keWeightedSpectra = GetKeWeightedLinealSpectra();

	//Open the file and write the value
	TFile* keWeightedLinealSpectraOutputFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraCellStudy.root","RECREATE");
	keWeightedLinealSpectraOutputFile->WriteObject(&keWeightedSpectra, "Lineal_energy_library");

	//Close the file
	keWeightedLinealSpectraOutputFile->Close();
}