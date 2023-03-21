#include "BWF_Plotter.h"
#include "TLine.h"

void GeneralizedBWFMultigraphPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BWF_Fitting_Results results, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<double>* const LETds = &survivalParams.LETd;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nAlphaparams = results.alphaFunc.GetNumFittingParams();
	const int nParams = (results.alphaFunc.GetNumFittingParams()+results.betaFunc.GetNumFittingParams());
	//const double beta = fitFuncParams[nParams];

	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	for(double LETd:*LETds) //First we iterate over each LET
	{
		c->cd(i+1);

		TGraph* gr = new TGraph();

		double alphaPredicted = results.alphaFunc.GetValue(fitFuncParams, LETd);
		double betaPredicted = results.betaFunc.GetValue(&(fitFuncParams[nAlphaparams]), LETd);

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(j))+((betaPredicted)*(j*j));
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

void AlphaBetaMultigraphResiduals(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BWF_Fitting_Results results,  double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nAlphaparams = results.alphaFunc.GetNumFittingParams();
	const int nParams = (results.alphaFunc.GetNumFittingParams()+results.betaFunc.GetNumFittingParams());
	int i = 0;

	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{

		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		TGraph* gr = new TGraph();

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaPredicted = 0;
		double betaPredicted = 0;

		double largestVal = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryAlphaVal = results.alphaFunc.GetValue(results.alphaFunc.GetFittingParams(),center);
			alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = results.betaFunc.GetValue(results.betaFunc.GetFittingParams(),center);
			betaPredicted += ryBetaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

		}

	   //For that experiment iterate through every datapoint
	   for (int j = 0; j < survivalParams.dose[i].size(); ++j)
	   {
	   		double calculatedSF = (alphaPredicted*(survivalParams.dose[i][j]))+(betaPredicted*(survivalParams.dose[i][j]*survivalParams.dose[i][j]));
	   		calculatedSF = std::exp(-calculatedSF);

	   		//Add the difference point to the graph
	   		double dataPoint = calculatedSF-survivalParams.survivingFraction[i][j];
	   		if (std::fabs(dataPoint) > largestVal) { largestVal = std::fabs(dataPoint); }
       		gr->AddPoint(survivalParams.dose[i][j],calculatedSF-survivalParams.survivingFraction[i][j]);
	   }


		//Draw
		gr->Draw("AP");

		//Set axes
		gr->SetTitle("");
		//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
		gr->GetYaxis()->SetTitle("P(SF)-SF");
		gr->GetXaxis()->SetTitle("Dose [Gy]");
		gr->GetXaxis()->CenterTitle(true);
		gr->GetYaxis()->CenterTitle(true);
		gr->GetXaxis()->SetTitleFont(42);
		gr->GetYaxis()->SetTitleFont(42);
		gr->GetXaxis()->SetTitleSize(0.052);
		gr->GetYaxis()->SetTitleSize(0.058);
		gr->GetXaxis()->SetTitleOffset(0.85);
		gr->GetYaxis()->SetTitleOffset(0.85);
		// gr->GetYaxis()->SetRangeUser(-0.5,0.5);
		gr->GetXaxis()->SetLimits(0,5.5);
		
		//Make the limit symmetric about Y
		double minAbs = std::fabs(c->GetUymin());
		double maxAbs = std::fabs(c->GetUymax());
		if (minAbs > maxAbs) {maxAbs = minAbs;} //See if the negative or positive side is larger
		gr->GetYaxis()->SetRangeUser(-largestVal*1.05,largestVal*1.05); //set the new y Limits

		//Set histogram fill and line settings
		gr->SetLineWidth(5);

		gr->SetMarkerColor(kRed);
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

		TLine* l = new TLine(0,0,5.5,0);
		l->SetLineStyle(9);
		l->SetLineWidth(4);
		l->Draw();

		++i;
	}
}

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
	const double beta = fitFuncParams[nParams];

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
		double betaPredicted = 0;

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

void MultigraphSurvivalFitPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, double* alphas_and_betas, double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;

	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	//Now let's iterate through our lineal energy library, plot, and do all that nice stuff
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
		c->cd(i+1);

		TGraph* gr = new TGraph();

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphas_and_betas[i*2]*(j))+((alphas_and_betas[(i*2)+1])*(j*j));
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

			// if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; }
			// legend->Draw();

			++i;
	}
}

void SurvivalDataMultigraphResiduals(TCanvas* c, TLegend* legend, CellStudyBWFFittingParameters& survivalParams, double* alphas_and_betas)
{
	//Iterate through each of the dose levels
	for (int i = 0; i < survivalParams.dose.size(); ++i)
	{
		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		TGraph* gr = new TGraph();
		double largestVal = 0;

	   //At that dose level iterate through every datapoint
	   for (int j = 0; j < survivalParams.dose[i].size(); ++j)
	   {
	   		double calculatedSF = (alphas_and_betas[(i*2)]*(survivalParams.dose[i][j]))+(alphas_and_betas[(i*2)+1]*(survivalParams.dose[i][j]*survivalParams.dose[i][j]));
	   		calculatedSF = std::exp(-calculatedSF);

	   		//Add the difference point to the graph
	   		double dataPoint = calculatedSF-survivalParams.survivingFraction[i][j];
	   		if (std::fabs(dataPoint) > largestVal) { largestVal = std::fabs(dataPoint); }
       		gr->AddPoint(survivalParams.dose[i][j],calculatedSF-survivalParams.survivingFraction[i][j]);
	   }


		//Draw
		gr->Draw("AP");

		//Set axes
		gr->SetTitle("");
		//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
		gr->GetYaxis()->SetTitle("P(SF)-SF");
		gr->GetXaxis()->SetTitle("Dose [Gy]");
		gr->GetXaxis()->CenterTitle(true);
		gr->GetYaxis()->CenterTitle(true);
		gr->GetXaxis()->SetTitleFont(42);
		gr->GetYaxis()->SetTitleFont(42);
		gr->GetXaxis()->SetTitleSize(0.052);
		gr->GetYaxis()->SetTitleSize(0.058);
		gr->GetXaxis()->SetTitleOffset(0.85);
		gr->GetYaxis()->SetTitleOffset(0.85);
		// gr->GetYaxis()->SetRangeUser(-0.5,0.5);
		gr->GetXaxis()->SetLimits(0,5.5);
		
		//Make the limit symmetric about Y
		double minAbs = std::fabs(c->GetUymin());
		double maxAbs = std::fabs(c->GetUymax());
		if (minAbs > maxAbs) {maxAbs = minAbs;} //See if the negative or positive side is larger
		gr->GetYaxis()->SetRangeUser(-largestVal*1.05,largestVal*1.05); //set the new y Limits

		//Set histogram fill and line settings
		gr->SetLineWidth(5);

		gr->SetMarkerColor(kRed);
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

		TLine* l = new TLine(0,0,5.5,0);
		l->SetLineStyle(9);
		l->SetLineWidth(4);
		l->Draw();
	}
}

void GeneralizedLETMultigraphPlotterAlphaBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<double>* const LETds = &survivalParams.LETd;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nAlphaparams = alphaFittingFunction.GetNumFittingParams();
	const int nParams = (alphaFittingFunction.GetNumFittingParams()+betaFittingFunction.GetNumFittingParams());
	//const double beta = fitFuncParams[nParams];

	bool legendAdded = false; 
	if (legendName == "") { legendAdded = true; }
	int i = 0;

	for(double LETd:*LETds) //First we iterate over each LET
	{
		c->cd(i+1);

		TGraph* gr = new TGraph();

		double alphaPredicted = alphaFittingFunction.GetValue(fitFuncParams, LETd);
		double betaPredicted = betaFittingFunction.GetValue(&(fitFuncParams[nAlphaparams]), LETd);

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(j))+((betaPredicted)*(j*j));
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

void GeneralizedBWFMultigraphPlotterAlphaBeta(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nAlphaparams = alphaFittingFunction.GetNumFittingParams();
	const int nParams = (alphaFittingFunction.GetNumFittingParams()+betaFittingFunction.GetNumFittingParams());
	//const double beta = fitFuncParams[nParams];

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
		double betaPredicted = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryAlphaVal = alphaFittingFunction.GetValue(fitFuncParams,center);
			alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = betaFittingFunction.GetValue(&fitFuncParams[nAlphaparams],center);
			betaPredicted += ryBetaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

		}

		for (double j = minDose; j < maxDose; j+=0.01) //Iterate through every dose at a given lineal energy
		{
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(j))+((betaPredicted)*(j*j));
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

void AlphaBetaMultigraphResiduals(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, const CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction alphaFittingFunction, BiologicalWeightingFunction betaFittingFunction, double* fitFuncParams,  double minDose, double maxDose)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	const int nAlphaparams = alphaFittingFunction.GetNumFittingParams();
	const int nParams = (alphaFittingFunction.GetNumFittingParams()+betaFittingFunction.GetNumFittingParams());
	int i = 0;


	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{

		c->cd(i+1);
		c->SetFillStyle(4000);
		c->SetFrameFillStyle(4000);

		TGraph* gr = new TGraph();

		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaPredicted = 0;
		double betaPredicted = 0;

		double largestVal = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryAlphaVal = alphaFittingFunction.GetValue(fitFuncParams,center);
			alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = betaFittingFunction.GetValue(&fitFuncParams[nAlphaparams],center);
			betaPredicted += ryBetaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

		}

	   //For that experiment iterate through every datapoint
	   for (int j = 0; j < survivalParams.dose[i].size(); ++j)
	   {
	   		double calculatedSF = (alphaPredicted*(survivalParams.dose[i][j]))+(betaPredicted*(survivalParams.dose[i][j]*survivalParams.dose[i][j]));
	   		calculatedSF = std::exp(-calculatedSF);

	   		//Add the difference point to the graph
	   		double dataPoint = calculatedSF-survivalParams.survivingFraction[i][j];
	   		if (std::fabs(dataPoint) > largestVal) { largestVal = std::fabs(dataPoint); }
       		gr->AddPoint(survivalParams.dose[i][j],calculatedSF-survivalParams.survivingFraction[i][j]);
	   }


		//Draw
		gr->Draw("AP");

		//Set axes
		gr->SetTitle("");
		//gr->SetTitleSize(0.03,"t"); //this doesn't do anything
		gr->GetYaxis()->SetTitle("P(SF)-SF");
		gr->GetXaxis()->SetTitle("Dose [Gy]");
		gr->GetXaxis()->CenterTitle(true);
		gr->GetYaxis()->CenterTitle(true);
		gr->GetXaxis()->SetTitleFont(42);
		gr->GetYaxis()->SetTitleFont(42);
		gr->GetXaxis()->SetTitleSize(0.052);
		gr->GetYaxis()->SetTitleSize(0.058);
		gr->GetXaxis()->SetTitleOffset(0.85);
		gr->GetYaxis()->SetTitleOffset(0.85);
		// gr->GetYaxis()->SetRangeUser(-0.5,0.5);
		gr->GetXaxis()->SetLimits(0,5.5);
		
		//Make the limit symmetric about Y
		double minAbs = std::fabs(c->GetUymin());
		double maxAbs = std::fabs(c->GetUymax());
		if (minAbs > maxAbs) {maxAbs = minAbs;} //See if the negative or positive side is larger
		gr->GetYaxis()->SetRangeUser(-largestVal*1.05,largestVal*1.05); //set the new y Limits

		//Set histogram fill and line settings
		gr->SetLineWidth(5);

		gr->SetMarkerColor(kRed);
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

		TLine* l = new TLine(0,0,5.5,0);
		l->SetLineStyle(9);
		l->SetLineWidth(4);
		l->Draw();

		++i;
	}
}


TGraph* BWFFunctionPlotter(TCanvas* c, TLegend* legend, const TAttLine& lineAttributes, std::string legendName, BiologicalWeightingFunction fittingFunction, double* fitFuncParams, std::string options, double minLineal, double maxLineal)
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
	gr->GetXaxis()->SetLimits(minLineal,maxLineal);

	//No marker
	gr->SetMarkerSize(0);

	if (!legendAdded) { legend->AddEntry(gr, (TString)legendName, "L"); legendAdded = true; legend->Draw();}
	

	return gr;

}

