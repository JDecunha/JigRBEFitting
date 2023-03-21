#include "BWF_Fitting_Results.h"

void BWF_Fitting_Results::PrintSummary()
{
	std::cout << summary.BriefReport() << std::endl;
}

void BWF_Fitting_Results::PrintBasic()
{
	std::cout << "Alpha fit params: "; alphaFunc.PrintFitParams();
	std::cout << std::endl << "Beta fit params: "; betaFunc.PrintFitParams();
	std::cout << std::endl << "RMSE: " << std::sqrt(2*summary.final_cost) << std::endl;

	std::cout << std::endl;
}