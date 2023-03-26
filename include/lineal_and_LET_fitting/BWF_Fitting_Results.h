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
		void PrintAIC(int const& numObservations);
		void PrintRMSEMinusPenaltyFunction(double penaltyWeight, int const& numObservations, double lower = 0.1, double upper = 250.);
};