#include "Ceres_BWF_Fitter.h"
#include "Ceres_LET_Fitter.h"



BWF_Fitting_Results Ceres_LET_Fitter::Fit()
{
	GeneralizedLETFitting();

	return results;
}

void Ceres_LET_Fitter::SetupResidualBlocks()
{
	int l = 0; //To keep track of jig position we are at 

	for(const std::pair<std::string,TH1D>& spectrumPair:_survivalParams.dySpectra) //First we iterate over each of the lineal energy spectra
	{
		TH1D const& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object

		for (int j = 0; j < _survivalParams.dose[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	
			ceres::DynamicNumericDiffCostFunction<Generalized_LET_Residual>* costFunc = new ceres::DynamicNumericDiffCostFunction<Generalized_LET_Residual>
				(new Generalized_LET_Residual(_survivalParams.dose[l][j], _survivalParams.survivingFraction[l][j], _survivalParams.LETd[l], _alphaFitFunc, _betaFitFunc));

			//Since we constructed a dynamic cost function, we have to tell it the number of parameters
			costFunc->AddParameterBlock(_alphaFitFunc.GetNumFittingParams());
			costFunc->AddParameterBlock(_betaFitFunc.GetNumFittingParams());
			costFunc->SetNumResiduals(1); //Number of residuals is just 1, because 1 SF is calculated at a time

			//Add a residual block with the alpha and beta BWF Fitting params
			_problem.AddResidualBlock(costFunc, nullptr, _alphaFitFunc.GetFittingParams(), _betaFitFunc.GetFittingParams()); 
		}

		++l; //iterate jig positions
	}

	//If positive constrained add the penalty term
	if (_positiveConstrained)
	{
		ceres::DynamicNumericDiffCostFunction<LET_Function_Negative_Penalty>* costFunc = new ceres::DynamicNumericDiffCostFunction<LET_Function_Negative_Penalty>
			(new LET_Function_Negative_Penalty(_alphaFitFunc, _betaFitFunc, _penaltyWeight));

		//Since we constructed a dynamic cost function, we have to tell it the number of parameters
		costFunc->AddParameterBlock(_alphaFitFunc.GetNumFittingParams());
		costFunc->AddParameterBlock(_betaFitFunc.GetNumFittingParams());
		costFunc->SetNumResiduals(1); //In this case, technically there are 2 residuals. But I combine them internally in the cost function

		//Add a residual block with the alpha and beta BWF Fitting params
		_problem.AddResidualBlock(costFunc, nullptr, _alphaFitFunc.GetFittingParams(), _betaFitFunc.GetFittingParams()); 
	}

}

void Ceres_LET_Fitter::GeneralizedLETFitting()
{	
	//If user has specified options, use them, if not use defaults
	ceres::Solver::Options options;
	if (_optionsSet) {options = _options;}
	else 
	{	
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = false;

		//Add the callback
		Ceres_BWF_Fitting_Callback callback(&_alphaFitFunc, &_betaFitFunc); //This callback is borrowed from the BWF fitter (it's a dependancy)
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

void Ceres_LET_Fitter::SetAlphaParameterLowerConstraint(int index, double constraint)
{
	_problem.SetParameterLowerBound(_alphaFitFunc.GetFittingParams(),index, constraint);
}

void Ceres_LET_Fitter::SetAlphaParameterUpperConstraint(int index, double constraint)
{
	_problem.SetParameterUpperBound(_alphaFitFunc.GetFittingParams(),index, constraint);
}

void Ceres_LET_Fitter::SetBetaParameterLowerConstraint(int index, double constraint)
{
	_problem.SetParameterLowerBound(_betaFitFunc.GetFittingParams(),index, constraint);
}

void Ceres_LET_Fitter::CheckFunctionNegativity(BiologicalWeightingFunction const& BWF, double lower, double upper)
{
        //Want to calculate the area under the function
        //And also calculate how negative it is
        double totalIntegralArea = 0.;
        double areaUnderFunction = 0.;
        double maxNegativeUnder = 0.;

        for(double center = lower; center <= upper; center+=0.1) //Iterate over the lineal energy spectrum
        {
          //Accumulate the predicted value of alpha as we integrate the function
          double ryVal = BWF.GetValue(center);
          if (ryVal > 0)
          {
            totalIntegralArea += ryVal*0.1;
          }
          else
          {
            areaUnderFunction += ryVal*0.1;
            totalIntegralArea += std::abs(ryVal)*0.1;
            if (ryVal < maxNegativeUnder) {maxNegativeUnder = ryVal;}
          }
        }

        if (maxNegativeUnder == 0.) { std::cout << "function is not negative." << std::endl; }
        else
        {
          std::cout << "Area under Y-axis: " << std::abs(areaUnderFunction) << std::endl;
          std::cout << "Fractional negative area: " << std::abs(areaUnderFunction*100.)/totalIntegralArea << "%" << std::endl;
          std::cout << "Maximum negative value: " << maxNegativeUnder << std::endl;
        }
}