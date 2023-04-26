#pragma once

#include "ceres/ceres.h"
#include "BiologicalWeightingFunction.h"
#include "NonLinearLeastSquaresUtilities.h"

class Survival_Fitting_Results
{
	public:
		ceres::Solver::Summary summary;
		std::vector<BiologicalWeightingFunction> alphabetaParams;

		void PrintSummary();
		void PrintBasic();
};