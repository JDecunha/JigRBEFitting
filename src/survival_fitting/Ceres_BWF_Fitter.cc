#include "Ceres_BWF_Fitter.h"

void Ceres_BWF_Fitter::Fit()
{
	if(!_paramsSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting survival params. Failure." << std::endl; return; }
	if(!_alphaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  alpha biological weighting function to fit. Failure." << std::endl; return; }
	if(!_betaFunctionSet) {std::cout << "Attempt to use Ceres_BWF_Fitter without setting  beta biological weighting function to fit. Failure." << std::endl; return; }

	GeneralizedBWFFitting();
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
	SetupResidualBlocks();

	// Run the solver!
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &_problem, &summary);

	std::cout << summary.BriefReport() << "\n";
}