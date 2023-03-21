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

void SurvivalDataMultigraph(TCanvas* c, TLegend* legend, CellStudySurvivalParameters survivalParams)
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