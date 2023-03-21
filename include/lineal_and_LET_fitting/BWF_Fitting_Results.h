#pragma once

#include "ceres/ceres.h"
#include "BiologicalWeightingFunction.h"

class BWF_Fitting_Results
{
	public:
		BiologicalWeightingFunction alphaFunc;
		BiologicalWeightingFunction betaFunc;
		ceres::Solver::Summary summary;

		void PrintSummary();
		void PrintBasic();
};