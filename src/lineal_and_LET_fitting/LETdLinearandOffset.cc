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
#include "TMultiGraph.h"
#include "TLegend.h"

//GNU Scientific Library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

//This project
#include "LETdLinearandOffset.h"
#include "NonLinearLeastSquaresUtilities.h"

int LinearSurvivalLETd(const gsl_vector* x, void* data, gsl_vector* f)
{
	//Load in the dySpectra and underlying ground truth alpha values
	const std::vector<double>* const LETds = &(((CellStudySurvivalParametersLETd*)data)->LETd);
	const std::vector<std::vector<double>>* const sf = &(((CellStudySurvivalParametersLETd*)data)->survivingFraction);
	const std::vector<std::vector<double>>* const dose = &(((CellStudySurvivalParametersLETd*)data)->dose);
	const double* beta = &(((CellStudySurvivalParametersLETd*)data)->beta);

	//Pull the current fitting parameters. For this function r(y) = c1*y+c0
	double c0 = gsl_vector_get(x,0);
	double c1 = gsl_vector_get(x,1);

	//I know, this is confusing because I have 4 different iterators
	// i for the overall # of iterations, k for the histogram bins, j for the dose level, l for the LET
	int i = 0; //To keep track of which LET and Dose we are on


	for(int l = 0; l < (*LETds).size(); ++l)
	{
		double alphaPredicted = 0;

		for (int j = 0; j < (*dose)[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	
			double survivalPredicted = 0; 

			alphaPredicted = (c1*(*LETds)[i])+c0;

			//calculate the exponent of survival predicted
			survivalPredicted = alphaPredicted*(((*dose)[l])[j])+(*beta)*(((*dose)[l])[j])*(((*dose)[l])[j]);
			//take e^exponent
			survivalPredicted = std::exp(-survivalPredicted);

			//Add entry to cost function vector with value Y_Predicted - Y_GroundTruth
			gsl_vector_set(f,i,survivalPredicted-((*sf)[l])[j]); 
			++i;
		}
	}
	
	return GSL_SUCCESS;
}

TCanvas* PlotSurvivalFromLinearLETModel(TCanvas* c, TLegend* legend, const std::vector<std::vector<double>> doselist, const std::vector<double> LETs, const double c0, const double c1, const double beta)
{
	bool legendAdded = false;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for (int i = 0; i < doselist.size(); ++i)
	{
			c->cd(i+1);

			TGraph* gr = new TGraph();

		   for (int j = 0; j < doselist[i].size(); ++j) //Iterate through every survival experiment at a given LET
		   {
		   		double alpha = (LETs[i]*c1)+c0;
		   		double calculatedSF = (alpha*(doselist[i][j]))+(beta*(doselist[i][j]*doselist[i][j]));
		   		calculatedSF = std::exp(-calculatedSF);
	       		gr->AddPoint(doselist[i][j],calculatedSF);
		   }
			//Draw
			gr->Draw("PL");

			//Set histogram fill and line settings
			gr->SetLineColor(kBlack);
			gr->SetLineWidth(0);
			gr->SetLineStyle(1);

			gr->SetMarkerColor(kBlue+2);
			gr->SetMarkerSize(5);
			gr->SetMarkerStyle(22);

			if (!legendAdded) { legend->AddEntry(gr,"#alpha = LET_{d} #upoint c_{1} + c_{0}","P"); legendAdded = true; }
			legend->Draw();
	}

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/multigraphSurvivingFractionWithLETdmodel.jpg";
	c->SaveAs((TString)outputName);

	return c;
}