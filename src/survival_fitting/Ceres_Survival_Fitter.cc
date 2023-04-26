#include "Ceres_Survival_Fitter.h"


Survival_Fitting_Results Ceres_Survival_Fitter::Fit()
{
	GeneralizedSurvivalFitting();

	return results;
}

void Ceres_Survival_Fitter::SetupResidualBlocks()
{
	int l = 0; //To keep track of jig position we are at 

	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	for(const std::pair<std::string,TH1D>& spectrumPair:_survivalParams.dySpectra) //First we iterate over each of the lineal energy spectra
	{
		results.alphabetaParams.push_back(FixedBWF); results.alphabetaParams.push_back(FixedBWF); 

		for (int j = 0; j < _survivalParams.dose[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	
			ceres::DynamicNumericDiffCostFunction<Generalized_Survival_Residual>* costFunc = new ceres::DynamicNumericDiffCostFunction<Generalized_Survival_Residual>
				(new Generalized_Survival_Residual(_survivalParams.dose[l][j], _survivalParams.survivingFraction[l][j]));

			//Add a parameter block for alpha and beta
			costFunc->AddParameterBlock(1);
			costFunc->AddParameterBlock(1);
			costFunc->SetNumResiduals(1); //Number of residuals is just 1, because 1 SF is calculated at a time

			//Add a residual block with the alpha and beta BWF Fitting params
			_problem.AddResidualBlock(costFunc, nullptr, results.alphabetaParams[l*2].GetFittingParams(), results.alphabetaParams[(l*2)+1].GetFittingParams()); 
		}

		++l; //iterate jig positions
	}

}

void Ceres_Survival_Fitter::GeneralizedSurvivalFitting()
{	
	//If user has specified options, use them, if not use defaults
	ceres::Solver::Options options;
	if (_optionsSet) {options = _options;}
	else 
	{	
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = false;

		//Add the callback
		Ceres_Survival_Fitting_Callback callback;
		options.callbacks.push_back(&callback);
		options.update_state_every_iteration = true;
	}

	ceres::Solver::Summary summary;
	ceres::Solve(options, &_problem, &summary);

	results.summary = summary;

	_alreadyRun = true; //Set the flag to indicate the fitter has run
}
