#include "BWF_Fitting_Results.h"
#include <math.h>

void BWF_Fitting_Results::PrintSummary()
{
	std::cout << summary.BriefReport() << std::endl;
}

void BWF_Fitting_Results::PrintBasic()
{
	//Calc RMSE
	double RMSE = std::sqrt(2*summary.final_cost);

	//Print
	std::cout << "Alpha fit params: "; alphaFunc.PrintFitParams();
	std::cout << std::endl << "Beta fit params: "; betaFunc.PrintFitParams();
	std::cout << std::endl << "RMSE: " << RMSE << std::endl;
}

void BWF_Fitting_Results::PrintAIC(double penaltyWeight, int const& numObservations, double lower, double upper)
{
	//Calc RMSE
	double RMSE = std::sqrt(2*summary.final_cost);
	double penaltyTerm = 0.;

	//Now to calculate penalty term
	for(double center = lower; center <= upper; center+=0.1) //Iterate over the lineal energy spectrum
    {
	    //Accumulate the predicted value of alpha as we integrate the function
		double ryAlphaVal = alphaFunc.GetValue(center);
		double ryBetaVal = betaFunc.GetValue(center);

		//Evaluate penalty terms
		if(ryAlphaVal < 0)
		{
			penaltyTerm += -ryAlphaVal*penaltyWeight;
		}
		//Evaluate penalty terms
		if(ryBetaVal < 0)
		{
			penaltyTerm += -ryBetaVal*penaltyWeight;
		}
     }

     //The cost function takes each term * 1/2, so to subtract out the penalty I have to subtract half of it
     double CostMinusPenalty = summary.final_cost-(penaltyTerm/2.);
     double RMSE_Corrected = std::sqrt(2*CostMinusPenalty);

 	//Print
	// std::cout << "Alpha fit params: "; alphaFunc.PrintFitParams();
	// std::cout << std::endl << "Beta fit params: "; betaFunc.PrintFitParams();
	// std::cout << std::endl << "RMSE: " << RMSE << std::endl;
    // std::cout << "RMSE Corrected: " << RMSE_Corrected << std::endl;

	//Calc RMSE
	int nParams = summary.num_parameters; //std::cout << "Num params: " << nParams << std::endl;
	double AIC = (numObservations*std::log(RMSE_Corrected*RMSE_Corrected))+2*nParams;
	std::cout << AIC << std::endl;
}

//Let's change this to print the penalty function and full RMSE seperately

void BWF_Fitting_Results::PrintRMSEMinusPenaltyFunction(double penaltyWeight, int const& numObservations, double lower, double upper)
{
	//Calc RMSE
	double RMSE = std::sqrt(2*summary.final_cost);
	double penaltyTerm = 0.;

	//Now to calculate penalty term
	for(double center = lower; center <= upper; center+=0.1) //Iterate over the lineal energy spectrum
    {
	    //Accumulate the predicted value of alpha as we integrate the function
		double ryAlphaVal = alphaFunc.GetValue(center);
		double ryBetaVal = betaFunc.GetValue(center);

		//Evaluate penalty terms
		if(ryAlphaVal < 0)
		{
			penaltyTerm += -ryAlphaVal*penaltyWeight;
		}
		//Evaluate penalty terms
		if(ryBetaVal < 0)
		{
			penaltyTerm += -ryBetaVal*penaltyWeight;
		}
     }

     //The cost function takes each term * 1/2, so to subtract out the penalty I have to subtract half of it
     double CostMinusPenalty = summary.final_cost-(penaltyTerm/2.);
     double RMSE_Corrected = std::sqrt(2*CostMinusPenalty);

 	//Print
	std::cout << "Alpha fit params: "; alphaFunc.PrintFitParams();
	std::cout << std::endl << "Beta fit params: "; betaFunc.PrintFitParams();
	std::cout << std::endl << "RMSE: " << RMSE << std::endl;
    std::cout << "RMSE Corrected: " << RMSE_Corrected << std::endl;

	//Calc RMSE
	int nParams = summary.num_parameters; std::cout << "Num params: " << nParams << std::endl;
	double AIC = (numObservations*std::log(RMSE_Corrected*RMSE_Corrected))+2*nParams;
	std::cout << "AIC: " << AIC << std::endl;
}

void BWF_Fitting_Results::PrintRMSEBSurvivingFractionandAIC(const CellStudyBWFFittingParameters& survivalParams, int const& numObservations)
{
	const std::vector<std::pair<std::string,TH1D>> DySpectra = survivalParams.dySpectra;
	const std::vector<std::vector<double>> doselist = survivalParams.dose;
	double RMSE = 0.;
	int i = 0;

	//Iterate through each of the experiments levels
	for(const std::pair<std::string,TH1D>& spectrumPair:DySpectra) //First we iterate over each of the lineal energy spectra
	{
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
			double ryAlphaVal = this->alphaFunc.GetValue(center);
			alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = this->betaFunc.GetValue(center);
			betaPredicted += ryBetaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

		}

	   //For that experiment iterate through every datapoint
	   for (int j = 0; j < survivalParams.dose[i].size(); ++j)
	   {
	   		double calculatedSF = (alphaPredicted*(survivalParams.dose[i][j]))+(betaPredicted*(survivalParams.dose[i][j]*survivalParams.dose[i][j]));
	   		calculatedSF = std::exp(-calculatedSF);

	   		RMSE += ((survivalParams.survivingFraction[i][j] - calculatedSF)*(survivalParams.survivingFraction[i][j] - calculatedSF));
	   }
		++i;
	}

	RMSE = std::sqrt(RMSE);

	std::cout << std::endl << "RMSE: " << RMSE << std::endl;

	//Calc AIC
	int nParams = summary.num_parameters; std::cout << "Num params: " << nParams << std::endl;
	double AIC = (numObservations*std::log(RMSE*RMSE))+2*nParams;
	std::cout << "AIC: " << AIC << std::endl;


}