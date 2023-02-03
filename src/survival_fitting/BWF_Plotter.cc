#include "BWF_Plotter.h"

void GeneralizedBWFMultigraphPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudySurvivalParameters& survivalParams, BiologicalWeightingFunction fittingFunction, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const double beta = survivalParams.beta;
	const int nParams = fittingFunction.GetNumFittingParams();

	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);

		TGraph* gr = new TGraph();

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaPredicted = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryVal = fittingFunction.GetValue(fitFuncParams,center);
			alphaPredicted += ryVal*value*width; //value*width is d(y)*dy, and everything else is r(y)
		}

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(j))+((beta)*(j*j));
			//take e^exponent
			survivalPredicted = std::exp(-survivalPredicted);
			//Add point to the graph
			gr->AddPoint(j,survivalPredicted);
		}
			//Draw
			gr->Draw("L");

			//Set line settings
			gr->SetLineColor(lineAttributes.GetLineColor());
			gr->SetLineWidth(lineAttributes.GetLineWidth());
			gr->SetLineStyle(lineAttributes.GetLineStyle());

			//No marker
			gr->SetMarkerSize(0);

			if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; }
			legend->Draw();

			++i;
	}
}

void GeneralizedBWFMultigraphPlotterBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudySurvivalParameters& survivalParams, BiologicalWeightingFunction fittingFunction, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nParams = fittingFunction.GetNumFittingParams();
	const double beta = fitFuncParams[nParams-1];

	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);

		TGraph* gr = new TGraph();

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaPredicted = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryVal = fittingFunction.GetValue(fitFuncParams,center);
			alphaPredicted += ryVal*value*width; //value*width is d(y)*dy, and everything else is r(y)
		}

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(j))+((beta)*(j*j));
			//take e^exponent
			survivalPredicted = std::exp(-survivalPredicted);
			//Add point to the graph
			gr->AddPoint(j,survivalPredicted);
		}
			//Draw
			gr->Draw("L");

			//Set line settings
			gr->SetLineColor(lineAttributes.GetLineColor());
			gr->SetLineWidth(lineAttributes.GetLineWidth());
			gr->SetLineStyle(lineAttributes.GetLineStyle());

			//No marker
			gr->SetMarkerSize(0);

			if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; }
			legend->Draw();

			++i;
	}
}

void BWFFunctionPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, BiologicalWeightingFunction fittingFunction, double* fitFuncParams, std::string options, double minLineal, double maxLineal)
{
	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	TGraph* gr = new TGraph();

	for (double j = minLineal; j < maxLineal; j+=0.1) //Iterate through every dose at a given lineal energy
	{
		gr->AddPoint(j,fittingFunction.GetValue(fitFuncParams,j));
	}

	//Draw
	gr->Draw(options.c_str());

	//Set line settings
	gr->SetLineColor(lineAttributes.GetLineColor());
	gr->SetLineWidth(lineAttributes.GetLineWidth());
	gr->SetLineStyle(lineAttributes.GetLineStyle());

	gr->GetYaxis()->SetTitle("r(y)");
	gr->GetXaxis()->SetTitle("y [#frac{keV}{#mum}]");
	gr->GetXaxis()->CenterTitle(true);
	gr->GetYaxis()->CenterTitle(true);
	gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitleFont(42);
	gr->GetXaxis()->SetTitleSize(0.042);
	gr->GetYaxis()->SetTitleSize(0.048);
	gr->GetXaxis()->SetTitleOffset(0.); //Offset x axis so no overlap
	gr->GetYaxis()->SetTitleOffset(0.); //Offset x axis so no overlap
	gr->GetXaxis()->SetLimits(0.,250.);

	//No marker
	gr->SetMarkerSize(0);

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; }
	legend->Draw();

}

