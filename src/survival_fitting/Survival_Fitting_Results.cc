#include "Survival_Fitting_Results.h"
#include <math.h>

void Survival_Fitting_Results::PrintSummary()
{
	for (int i = 0; i < alphabetaParams.size(); i=i+2)
	{
		std::cout << "Alpha: " << alphabetaParams[i].GetValue(0) << " Beta: " << alphabetaParams[i+1].GetValue(0) << std::endl;
	}

	std::cout << summary.BriefReport() << std::endl;
}

void Survival_Fitting_Results::PrintBasic()
{
	
}