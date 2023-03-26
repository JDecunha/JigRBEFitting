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

void BWF_Fitting_Results::PrintAIC(int const& numObservations)
{
	//Calc RMSE
	double RMSE = std::sqrt(2*summary.final_cost);
	int nParams = summary.num_parameters; std::cout << "Num params: " << nParams << std::endl;
	double AIC = (numObservations*std::log(RMSE*RMSE))+2*nParams;
	std::cout << "AIC: " << AIC << std::endl;
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