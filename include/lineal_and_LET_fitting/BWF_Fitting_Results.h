#pragma once

#include "ceres/ceres.h"
#include "BiologicalWeightingFunction.h"
#include "NonLinearLeastSquaresUtilities.h"

class BWF_Fitting_Results
{
	public:
		BiologicalWeightingFunction alphaFunc;
		BiologicalWeightingFunction betaFunc;
		ceres::Solver::Summary summary;

		void PrintSummary();
		void PrintBasic();
		void PrintAIC(int const& numObservations);
		void PrintRMSEMinusPenaltyFunction(double penaltyWeight, int const& numObservations, double lower = 0.1, double upper = 120.);

		//This prints the RMSE as calculated from the surviving fraction alone, not the log(SF) terms.
		void PrintRMSEBSurvivingFractionandAIC(const CellStudyBWFFittingParameters& survivalParams, int const& numObservations);
};