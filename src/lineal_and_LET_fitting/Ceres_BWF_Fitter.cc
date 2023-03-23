#include "Ceres_BWF_Fitter.h"



BWF_Fitting_Results Ceres_BWF_Fitter::Fit()
{
	GeneralizedBWFFitting();

	return results;
}

void Ceres_BWF_Fitter::SetupResidualBlocks()
{
	int l = 0; //To keep track of jig position we are at 

	for(const std::pair<std::string,TH1D>& spectrumPair:_survivalParams.dySpectra) //First we iterate over each of the lineal energy spectra
	{
		TH1D const& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object

		for (int j = 0; j < _survivalParams.dose[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	

			//Construct the cost function with the Dose, SF, dy spectrum, and fitting functions
			ceres::DynamicNumericDiffCostFunction<Generalized_BWF_Residual>* costFunc = new ceres::DynamicNumericDiffCostFunction<Generalized_BWF_Residual>
				(new Generalized_BWF_Residual(_survivalParams.dose[l][j], _survivalParams.survivingFraction[l][j], dySpectrum, _alphaFitFunc, _betaFitFunc));

			//Since we constructed a dynamic cost function, we have to tell it the number of parameters
			costFunc->AddParameterBlock(_alphaFitFunc.GetNumFittingParams());
			costFunc->AddParameterBlock(_betaFitFunc.GetNumFittingParams());
			costFunc->SetNumResiduals(1); //Number of residuals is just 1, because 1 SF is calculated at a time

			//Add a residual block with the alpha and beta BWF Fitting params
			_problem.AddResidualBlock(costFunc, nullptr, _alphaFitFunc.GetFittingParams(), _betaFitFunc.GetFittingParams()); 
		}

		++l; //iterate jig positions
	}

}

void Ceres_BWF_Fitter::GeneralizedBWFFitting()
{	
	//If user has specified options, use them, if not use defaults
	ceres::Solver::Options options;
	if (_optionsSet) {options = _options;}
	else 
	{	
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = false;

		//Add the callback
		Ceres_BWF_Fitting_Callback callback(&_alphaFitFunc, &_betaFitFunc);
		options.callbacks.push_back(&callback);
		options.update_state_every_iteration = true;
	}

	ceres::Solver::Summary summary;
	ceres::Solve(options, &_problem, &summary);

	results.alphaFunc = _alphaFitFunc;
	results.betaFunc = _betaFitFunc;
	results.summary = summary;

	_alreadyRun = true; //Set the flag to indicate the fitter has run
}

void Ceres_BWF_Fitter::SetAlphaParameterLowerConstraint(int index, double constraint)
{
	_problem.SetParameterLowerBound(_alphaFitFunc.GetFittingParams(),index, constraint);
}

void Ceres_BWF_Fitter::SetAlphaParameterUpperConstraint(int index, double constraint)
{
	_problem.SetParameterUpperBound(_alphaFitFunc.GetFittingParams(),index, constraint);
}

void Ceres_BWF_Fitter::SetBetaParameterLowerConstraint(int index, double constraint)
{
	_problem.SetParameterLowerBound(_betaFitFunc.GetFittingParams(),index, constraint);
}