//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include "SurvivalPlotting.h"
#include "NonLinearLeastSquaresUtilities.h"
#include "BWF_Fitting_Results.h"

void PlotAlphaBeta(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetas, bool plotBeta)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlot;
	if (!plotBeta)
	{
		for (int i = 0; i < survivalParams.LETd.size(); ++i)
		{
			toPlot.push_back(alphaBetas[i*2]);
		}
	}
	else
	{
		for (int i = 0; i < survivalParams.LETd.size(); ++i)
		{
			toPlot.push_back(alphaBetas[(i*2)+1]);
		}
	}

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(survivalParams.LETd.size(),survivalParams.LETd.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	if (!plotBeta) { gr->GetYaxis()->SetTitle("#alpha"); }
	else { gr->GetYaxis()->SetTitle("#beta"); }
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0.9);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void PlotAlphaBetaFromBWF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, BWF_Fitting_Results results, bool plotBeta)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }


	BiologicalWeightingFunction weightFunc;
	if (!plotBeta)
	{
		weightFunc = results.alphaFunc;
	}
	else
	{
		weightFunc = results.betaFunc;
	}

	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	int i = 0;

	TGraph* gr = new TGraph();

	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double bwfValue = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryVal = weightFunc.GetValue(weightFunc.GetFittingParams(),center);
			bwfValue += ryVal*value*width; //value*width is d(y)*dy, and everything else is r(y)
		}
		gr->AddPoint(survivalParams.LETd[i],bwfValue);
		++i;
	}

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	if (!plotBeta) { gr->GetYaxis()->SetTitle("#alpha"); }
	else { gr->GetYaxis()->SetTitle("#beta"); }
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0.9);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudyBWFFittingParameters survivalParams)
{
	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (int i = 0; i < survivalParams.dose.size(); ++i)
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		TGraphErrors* gr = new TGraphErrors(survivalParams.dose[i].size(),&survivalParams.dose[i][0],&survivalParams.survivingFraction[i][0], nullptr, &survivalParams.survivingFractionUncertainty[i][0]);

		//Draw
		gr->Draw("AP");

		//Set axes
		gr->SetTitle("");
		//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
		gr->GetYaxis()->SetTitle("SF");
		gr->GetXaxis()->SetTitle("Dose [Gy]");
		gr->GetXaxis()->CenterTitle(true);
		gr->GetYaxis()->CenterTitle(true);
		gr->GetXaxis()->SetTitleFont(42);
		gr->GetYaxis()->SetTitleFont(42);
		gr->GetXaxis()->SetTitleSize(0.052);
		gr->GetYaxis()->SetTitleSize(0.058);
		gr->GetYaxis()->SetTitleOffset(0.9);
		gr->GetYaxis()->SetRangeUser(0.01,0.99);
		gr->GetXaxis()->SetLimits(0,5.5);
		gPad->SetLogy();

		//Set histogram fill and line settings
		gr->SetLineWidth(5);

		gr->SetMarkerColor(kBlack);
		gr->SetMarkerSize(4);
		gr->SetMarkerStyle(8);

		// Create an output string stream
		std::ostringstream outputstream;
		// Set Fixed -Point Notation
		outputstream << std::fixed;
		// Set precision to 2 digits
		outputstream << std::setprecision(1);
		//Add double to stream
		outputstream << survivalParams.LETd[i];
		std::string titlestring = outputstream.str() + " keV/#mum";
		TLatex *t = new TLatex(0.015, 0.935, (TString)titlestring);
		t->SetNDC(); //set position in coordinate system relative to canvas
		t->Draw();
	}
}

void PlotSurvivalFromAlphaBeta(TCanvas* c, const std::vector<std::vector<double>> doselist, const std::vector<double> alphas, const std::vector<double> betas)
{
	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (int i = 0; i < doselist.size(); ++i)
	{
			c->cd(i+1);

			TGraph* gr = new TGraph();

		   for (int j = 0; j < doselist[i].size(); ++j)
		   {
		   		double calculatedSF = (alphas[i]*(doselist[i][j]))+(betas[i]*(doselist[i][j]*doselist[i][j]));
		   		calculatedSF = std::exp(-calculatedSF);
	       		gr->AddPoint(doselist[i][j],calculatedSF);
		   }

			//Draw
			gr->Draw("PL");

			//Set histogram fill and line settings
			gr->SetLineColor(kBlack);
			gr->SetLineWidth(0);
			gr->SetLineStyle(1);

			gr->SetMarkerColor(kAzure+6);
			gr->SetMarkerSize(5);
			gr->SetMarkerStyle(22);
			

	}

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/multigraphSurvivingFractionWithModel.jpg";
	c->SaveAs((TString)outputName);

	delete c;
}

void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, double const* alphaBetasProton)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlot;

	for (int i = 0; i < survivalParams.LETd.size(); ++i)
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
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(survivalParams.LETd.size(),survivalParams.LETd.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0.7);
	

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());
	gr->SetLineColor(lineAttributes.GetLineColor());
	gr->SetLineWidth(lineAttributes.GetLineWidth());
	gr->SetLineStyle(lineAttributes.GetLineStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, BWF_Fitting_Results results)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	int i = 0;

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlot;
	std::vector<double> alphaBetasProton;

	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
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

	

	for (int i = 0; i < survivalParams.LETd.size(); ++i)
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
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(survivalParams.LETd.size(),survivalParams.LETd.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
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

void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, std::vector<BiologicalWeightingFunction> alphabetaParams)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	int i = 0;

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlot;
	std::vector<double> alphaBetasProton;

	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);
		
		alphaBetasProton.push_back(alphabetaParams[2*i].GetValue(1.));
		alphaBetasProton.push_back(alphabetaParams[(2*i)+1].GetValue(1.));

		++i;
	}

	

	for (int i = 0; i < survivalParams.LETd.size(); ++i)
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
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(survivalParams.LETd.size(),survivalParams.LETd.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
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

void PlotRBE10SFLET(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, BWF_Fitting_Results results)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	const std::vector<double> LETDs = survivalParams.LETd;
	int i = 0;

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlot;
	std::vector<double> alphaBetasProton;

	//Iterate through each of the experiments levels
	for(const double _LETd:LETDs) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

	    double alphaPredicted = results.alphaFunc.GetValue(_LETd);
	    double betaPredicted = results.betaFunc.GetValue(_LETd);
		
		alphaBetasProton.push_back(alphaPredicted);
		alphaBetasProton.push_back(betaPredicted);

		++i;
	}

	

	for (int i = 0; i < survivalParams.LETd.size(); ++i)
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
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(survivalParams.LETd.size(),survivalParams.LETd.data(),toPlot.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0);
	gr->GetXaxis()->SetTitleOffset(1.15);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "P"); legendAdded = true; legend->Draw(); }
}

void PlotRBE10SFMcNamara(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes,  std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	//Calculate D_10% for cesium
	//quadratic formula
	double b = alphaBetasCesium[0];
	double a = alphaBetasCesium[1];
	double cval = -1;
	double D10Cesium = (-b+std::sqrt((b*b)-(4*a*cval)))/(2*a);

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlotx;std::vector<double> toPloty;

	//Steps:
	//Step 1.) Calculate RBE McNamara at all D_ps less than D_cesium for 0.1 SF
	//Step 2.) Iterate through the list of RBE McNamara and multiply D_p * RBE to give D_x
	//Step 3.) Find the D_x which is closest to D_x 0.1 SF and extract D_p from that.
	//Step 4.) Take RBE at that value

	double alphaBetaCesium = alphaBetasCesium[0]/alphaBetasCesium[1];
	double alphaBetaCesiumSquared = alphaBetaCesium*alphaBetaCesium;
	double p0 = 0.999064;
	double p1 = 0.35605;
	double p2 = 1.1012;
	double p3 = -0.0038703;

	for (double i = 0.8; i < 20; i += 1) //int i = 0; i < survivalParams.LETd.size(); ++i)
	{
		std::vector<double> ProtonDoses; 
		std::vector<double> RBEMcNamara; 
		std::vector<double> RBETimesProtonDose; 
		double LETd = i;//survivalParams.LETd[i];
		// std::cout << "LET: " << LETd << std::endl;
		// std::cout << "Alpha beta cesium ratio: " << alphaBetaCesium << std::endl;

		for (double j = 0.01; j <= D10Cesium; j+=0.01) //Calculate RBE McNamara at all the doses
		{
			double outFrontTerm = 1/(2*j);
			double term1 = alphaBetaCesiumSquared;
			double term2 = (4*j*alphaBetaCesium)*(p0+((p1/alphaBetaCesium)*LETd));
			double term3 = (4*j*j)*(p2-((p3*std::sqrt(alphaBetaCesium))*LETd))*(p2-((p3*std::sqrt(alphaBetaCesium))*LETd));
			
			double sqrtTerm = std::sqrt(term1+term2+term3);
			double RBE = outFrontTerm*(sqrtTerm-alphaBetaCesium);

			// std::cout << outFrontTerm << " " << term1 << " " << term2 << " " << term3 << std::endl;

			RBEMcNamara.push_back(RBE);
			ProtonDoses.push_back(j);
			RBETimesProtonDose.push_back(RBE*j);
			// std::cout << "Dp: " << j << " RBE: " << RBE << std::endl << std::endl;
		}

		double DiffDxProtonDx = 100;
		double DiffDxProtonDxLastIteration = 100;
		double RBE = 0.;

		for (int k  = 0; k < ProtonDoses.size(); ++k)
		{
			DiffDxProtonDxLastIteration = DiffDxProtonDx;
			DiffDxProtonDx = std::abs(RBETimesProtonDose[k]-D10Cesium);
			if (DiffDxProtonDx > DiffDxProtonDxLastIteration)
			{
				RBE = RBEMcNamara[k];
				// std::cout << RBE << std::endl;
				break;
			}
		}

		//Push back the ratio of doses
		toPlotx.push_back(i);
		toPloty.push_back(RBE);
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(toPlotx.size(),toPlotx.data(),toPloty.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0.65);
	gr->GetYaxis()->SetRangeUser(0.8,4); //set the new y Limits

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());
	gr->SetLineColor(lineAttributes.GetLineColor());
	gr->SetLineWidth(lineAttributes.GetLineWidth());
	gr->SetLineStyle(lineAttributes.GetLineStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; legend->Draw(); }
}

void PlotRBE10SFChenAndAhmad(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	int i = 0;

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlotx;
	std::vector<double> toPloty;
	std::vector<double> alphaBetasProton;
	

	for (double i = 0.8; i < 20; i += 1)
	{
		//quadratic formula
		double b = alphaBetasCesium[0];
		double a = alphaBetasCesium[1];
		double c = -1;
		double resultCesium = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Fitting params from Chen and Ahmad 2010
		double lambda1 = 0.0013;
		double lambda2 = 0.045;
		double LET2 = i*i;
		double top = 1-std::exp(-lambda1*LET2);
		double bottom = lambda2*i;
		b += (top/bottom);
		double resultProton = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Push back the ratio of doses
		toPloty.push_back(resultCesium/resultProton);
		toPlotx.push_back(i);
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(toPlotx.size(),toPlotx.data(),toPloty.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0);
	gr->GetXaxis()->SetTitleOffset(1.15);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());
	gr->SetLineColor(lineAttributes.GetLineColor());
	gr->SetLineWidth(lineAttributes.GetLineWidth());
	gr->SetLineStyle(lineAttributes.GetLineStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; legend->Draw(); }
}

void PlotRBE10SFWedenberg(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, const TAttLine& lineAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }

	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	int i = 0;

	//Assign either the alphas, or the betas to toPlot depending on the option
	std::vector<double> toPlotx; std::vector<double> toPloty;
	std::vector<double> alphaBetasProton;

	for (double i = 0.8; i < 20; i += 1)
	{
		//quadratic formula
		double b = alphaBetasCesium[0];
		double a = alphaBetasCesium[1];
		double c = -1;
		double resultCesium = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Fitting param from Wedenberg
		double q = 0.434;
		double alphaBetaRatio = b/a; //confusingly, b in the quadratic formula, is alpha
		double rightTerm = 1.+((q*i)/(alphaBetaRatio));
		b = rightTerm*alphaBetasCesium[0];
		double resultProton = (-b+std::sqrt((b*b)-(4*a*c)))/(2*a);

		//Push back the ratio of doses
		toPloty.push_back(resultCesium/resultProton);
		toPlotx.push_back(i);
	}
	

	//Constructor: Size and then two doubles
	TGraph* gr = new TGraph(toPlotx.size(),toPlotx.data(),toPloty.data());

	//Draw
	gr->Draw((TString)options);

	//Set axes
	gr->SetTitle("");
	//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
	gr->GetYaxis()->SetTitle("RBE(0.1 SF)"); 
	gr->GetXaxis()->SetTitle("LET_{d} [keV/#mum]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.052);
	gr->GetYaxis()->SetTitleSize(0.058);
	gr->GetYaxis()->SetTitleOffset(0);
	gr->GetXaxis()->SetTitleOffset(1.15);

	gr->SetMarkerColor(markerAttributes.GetMarkerColor());
	gr->SetMarkerSize(markerAttributes.GetMarkerSize());
	gr->SetMarkerStyle(markerAttributes.GetMarkerStyle());
	gr->SetLineColor(lineAttributes.GetLineColor());
	gr->SetLineWidth(lineAttributes.GetLineWidth());
	gr->SetLineStyle(lineAttributes.GetLineStyle());

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; legend->Draw(); }
}