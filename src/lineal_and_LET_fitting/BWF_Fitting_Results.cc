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