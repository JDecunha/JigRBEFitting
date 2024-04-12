//This project
#include "LinealSpectra.h"
#include "Utilities.h"

#include <iostream>
#include <fstream>

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

	std::vector<std::pair<std::string, TH1D>> monoenergeticSpectra; //List to hold our d(y) spectrum
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

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetKeWeightedFrequencyLinealSpectra(std::string targetSize)
{
	/*
	Retrieves the lineal energy spectra calculated at every location of the jig.
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	// std::string fyFolder = "/home/joseph/Documents/PHD_local/July_2022/proton_5umVoxel_DNA2_10kTracks";
	std::string fyFolder = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/5um_0to100MeV_April2023Library";


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
		TH1D outputLinealSpectrum = utils::GetFy(fyFolder, 0.1, targetSize);
		outputLinealSpectrum.Scale(0); //to eliminate the data and leave the bins

		//Loop over the KE spectrum and weight the f(y) spectra accordingly
		for (int i = 2; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			//Because of the way we re-defined the KE histogram, "GetBinLowEdge" actually corresponds to the middle
			TH1D tempLinealSpectrum = utils::GetNy(fyFolder, KESpectrum.GetBinLowEdge(i), targetSize); //Get the N(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Scale the spectrum so N(y) matches an effective count of 10 million
			long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, KESpectrum.GetBinLowEdge(i), targetSize);
			double scalingFactor = double(1e7)/double(effectiveNumTracks);
			tempLinealSpectrum.Scale(scalingFactor);

			//Add the weighted spectrum to the output fy spectrum
			outputLinealSpectrum.Add(&tempLinealSpectrum);
		}	
		
		//Convert N(y) to f(y)
		utils::PMF_to_FrequencyFunction(&outputLinealSpectrum);
		//utils::VerifyNormalization(outputLinealSpectrum); //Shows all spectra are normalized to 1 within 1e-7
		linealEnergyLibrary.push_back(std::make_pair<std::string,TH1D>(std::move(KEColumnName),std::move(outputLinealSpectrum)));
	}

	delete keSpectra; //we own the keSpectra pointer so we have to delete it

	return linealEnergyLibrary;
}

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetKeWeightedLinealSpectra(std::string targetSize)
{
	/*
	Retrieves the lineal energy spectra calculated at every location of the jig.
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Data/5um_0to100MeV_April2023Library";


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
			TH1D tempLinealSpectrum = utils::GetNy(fyFolder, KESpectrum.GetBinLowEdge(i), targetSize); //Get the N(y) spectrum at that energy

			//Scale the spectrum by the fluence
			tempLinealSpectrum.Scale(KESpectrum.GetBinContent(i));

			//Scale the spectrum so N(y) matches an effective count of 10 million
			long long effectiveNumTracks = utils::GetEffectiveNumberOfTracks(fyFolder, KESpectrum.GetBinLowEdge(i), targetSize);
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

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetKeWeightedLinealFrequencySpectraJuly2023(std::string targetSize)
{
	/*
	Retrieves the lineal energy spectra calculated at every location of the jig.
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/SuperTrack_01to100MeV_June2023_binspanning";

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<std::string,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/FadasProtonSpectra_July2023.root");
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

		//Loop over the KE spectrum and weight the f(y) spectra accordingly
		for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			//Because of the way we re-defined the KE histogram, "GetBinLowEdge" actually corresponds to the middle
			// std::cout << KESpectrum.GetBinCenter(i) << std::endl;
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

std::vector<std::pair<std::string,TH1D>> LinealSpectra::GetKeWeightedLinealSpectraJuly2023(std::string targetSize)
{
	/*
	Retrieves the lineal energy spectra calculated at every location of the jig.
	Returns a series of d(y) with jig location identifier.
	*/

	//Specify hardcoded path to our f(y) library and target size
	std::string fyFolder = "/home/joseph/Documents/PHD_local/SuperTrack_01to100MeV_June2023_binspanning";

	//Have to make a dictionary so that CERN ROOT can properly load and save STD library types from files
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");

	//Pull the data from our KE-spectrum file
	std::vector<std::pair<std::string,TH1F>>* keSpectra;
	TFile* spectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/FadasProtonSpectra_July2023.root");
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

		//Loop over the KE spectrum and weight the f(y) spectra accordingly
		for (int i = 1; i <= KESpectrum.GetNbinsX(); ++i) //start at 2 to skip underflow and 0.0 MeV bin
		{
			//Because of the way we re-defined the KE histogram, "GetBinLowEdge" actually corresponds to the middle
			// std::cout << KESpectrum.GetBinCenter(i) << std::endl;
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
			// utils::PMF_to_DoseFunction(spectra);
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

	std::vector<double> LETd{0.914, 1.16, 1.61, 1.81, 1.93, 2.34, 3.01, 5.08, 10.8, 15.2, 17.7, 19};

	int i = 1;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (const auto& linealSpectra:linealLibrary)
	{
			c->cd(i); ++i;//Change the window of the canvas we're drawing on (since it's a multigraph)

			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			// utils::PMF_to_DoseFunction(spectra);
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
			// Create an output string stream
			std::ostringstream outputstream;
			// Set Fixed -Point Notation
			outputstream << std::fixed;
			// Set precision to 2 digits
			outputstream << std::setprecision(1);
			//Add double to stream
			outputstream << LETd[i-2];
			std::string titlestring = outputstream.str() + " keV/#mum";
			TLatex *t = new TLatex(0.015, 0.935, (TString)titlestring);
			// TLatex *t = new TLatex(0.015, 0.935, (TString)(LETd[i-2]+ " keV/#mum"));
			t->SetNDC(); //set position in coordinate system relative to canvas
			t->Draw();
	}

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/linealSpectra/multigraphDyJigLinealJuly2023.jpg";
	c->SaveAs((TString)outputName);

	delete c;
}

void LinealSpectra::CompareGeant4toFluenceWeighting()
{
	gStyle->SetOptStat(0); //Don't print the stats window in the top right

	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(2040, 1640);
	c->SetWindowSize(2040, 1640);
	c->SetLeftMargin(0.15);
	c->SetBottomMargin(0.15);

	//Fada A
	// std::string path = "/home/joseph/Downloads/STV/proton_fada_A_1um_688576084.root";

	//Fada
	std::string path = "/home/joseph/Downloads/STV/proton_fada_G_1um_1501494966.root";
	
	//Get the KeWeighted lineal energy spectra
	auto keWeightedSpectra = LinealSpectra::GetKeWeightedLinealSpectra("1e3");

	TFile f = TFile((TString)path);
	//TFile f2 = TFile((TString)path2);

	TH1D* h;
	TH1D* h2 = new TH1D(std::get<1>(keWeightedSpectra[6]));
	h = (TH1D*)f.Get("f(y)");
	// h2 = (TH1D*)f2.Get("Lineal energy histogram");
	//Print_TH1(h); //To manually print the bin values
	Utilities::PMF_to_DoseFunction(h);
	Utilities::Prepare_for_Semilog(h);
	Utilities::PMF_to_DoseFunction(h2);
	Utilities::Prepare_for_Semilog(h2);

	gStyle->SetTitleFont(42,"t");

	h->Draw("HIST");
	h2->Draw("same");
	
	//Set range
	//h2->GetYaxis()->SetRangeUser(0,1);
	//h->GetYaxis()->SetRangeUser(0,1);

	//Set titles
	//std::string title = std::to_string(int(energy));
	//title += " MeV";
	h->SetTitle("");
	//h2->SetTitle("");
	h->SetTitleSize(0.03,"t"); //this doesn't do anything
	h->GetYaxis()->SetTitle("y #upoint d(y)");
	h->GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");

	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);
	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetTitleSize(0.042);
	h->GetYaxis()->SetTitleSize(0.048);
	h->GetXaxis()->SetTitleOffset(1.50); //Offset x axis so no overlap

	//Center
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);

	//Offset x axis so no overlap
	h->GetXaxis()->SetTitleOffset(1.2);
	h->GetYaxis()->SetTitleOffset(1.2);

	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);


	//h->SetFillColorAlpha(kRed, 0.4);
	//h->SetFillStyle(3001);
	h->SetLineColor(kBlack);
	h->SetLineWidth(7);

	//h2->SetFillColorAlpha(kAzure+3, 0.5);
	//h->SetFillStyle(3001);
	h2->SetFillColorAlpha(kAzure+3, 0.5);
	h2->SetLineColor(kBlack);
	h2->SetLineWidth(0);

	gPad->SetLogx();

	//auto legend = new TLegend(0.18,0.73,0.18+0.26,0.73+0.15);//x start, y start, x end, yend
	double xstart = 0.71;
	auto legend = new TLegend(xstart,0.73,xstart+0.17,0.73+0.15);//x start, y start, x end, yend
	//legend->SetHeader("","C"); // option "C" allows to center the header
	legend->AddEntry(h,"Geant4","L");
	legend->AddEntry(h2,"SuperTrack","f");
	//legend->AddEntry(g5,"10 #mum","L");
	//legend->AddEntry(g6,"1","L");
	//legend->AddEntry(g5,"1 Million Tracks","P");
	legend->SetTextSize(0.030);
	legend->Draw();

	c->Update();
	c->Print("/home/joseph/Dropbox/Documents/Work/Projects/MDA_Microdosimetry/Images/SuperTrackValidation/FadaG2.jpg"); 
}

void LinealSpectra::SaveKeWeightedLinealSpectra()
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& value : targetSizes)
	{
		auto keWeightedSpectra = GetKeWeightedLinealSpectraJuly2023(value);

		//Open the file and write the value
		TFile* keWeightedLinealSpectraOutputFile = TFile::Open((TString)("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraJuly2023CellStudy_"+value+"nm.root"),"RECREATE");
		keWeightedLinealSpectraOutputFile->WriteObject(&keWeightedSpectra, "Lineal_energy_library");

		//Close the file
		keWeightedLinealSpectraOutputFile->Close();
	}
}

void LinealSpectra::SaveKeWeightedLinealSpectraCSV()
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& value : targetSizes)
	{
		auto keWeightedSpectra = GetKeWeightedLinealSpectraJuly2023(value);

		for (const auto& linealSpectra:keWeightedSpectra)
		{
			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			auto columnName = std::get<0>(linealSpectra);

			std::ofstream myfile;
			myfile.open ("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/dyCSV/dySpectrum_"+columnName+".csv");
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

void LinealSpectra::SaveKeWeightedLinealFrequencySpectraCSV()
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& value : targetSizes)
	{
		auto keWeightedSpectra = GetKeWeightedLinealFrequencySpectraJuly2023(value);

		for (const auto& linealSpectra:keWeightedSpectra)
		{
			//Get the spectra and its ID
			TH1D* spectra = new TH1D(std::get<1>(linealSpectra));
			auto columnName = std::get<0>(linealSpectra);

			std::ofstream myfile;
			myfile.open ("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/fyCSV/fySpectrum_"+columnName+".csv");
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

void LinealSpectra::SaveKeWeightedFrequencyLinealSpectra()
{
	//Add entry to the dictionary so ROOT can read and save files with the given std library type
	gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");

	//Get the KeWeighted lineal energy spectra
	// std::vector<std::string> targetSizes{"10","50","100","200","300","400","500","600","700","800","900","1e3"};
	std::vector<std::string> targetSizes{"1e3"};
	for (const auto& value : targetSizes)
	{
		auto keWeightedSpectra = GetKeWeightedFrequencyLinealSpectra(value);

		//Open the file and write the value
		TFile* keWeightedLinealSpectraOutputFile = TFile::Open((TString)("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealFYSpectraApril2023CellStudy_"+value+"nm.root"),"RECREATE");
		keWeightedLinealSpectraOutputFile->WriteObject(&keWeightedSpectra, "Lineal_energy_library");

		//Close the file
		keWeightedLinealSpectraOutputFile->Close();
	}
}