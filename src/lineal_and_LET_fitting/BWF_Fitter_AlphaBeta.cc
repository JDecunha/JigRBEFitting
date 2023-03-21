#include "BWF_Fitter_AlphaBeta.h"

double* BWF_Fitter_AlphaBeta::Fit(double* initialGuess, bool weightedFit)
{
	if(!_paramsSet) {std::cout << "Attempt to use BWF_Fitter_AlphaBeta without setting survival params. Failure." << std::endl; return nullptr; }
	if(!_alphaFunctionSet) {std::cout << "Attempt to use BWF_Fitter_AlphaBeta without setting  alpha biological weighting function to fit. Failure." << std::endl; return nullptr; }
	if(!_betaFunctionSet) {std::cout << "Attempt to use BWF_Fitter_AlphaBeta without setting  beta biological weighting function to fit. Failure." << std::endl; return nullptr; }

	return GeneralizedBWFFitting(_survivalParams, _alphaFitFunc, _betaFitFunc, initialGuess, weightedFit);
}

int BWF_Fitter_AlphaBeta::GeneralizedBWF(const gsl_vector* x, void* data, gsl_vector* f)
{
	//Load in the dySpectra and underlying ground truth alpha values
	const std::vector<std::pair<std::string,TH1D>>* const pDySpectra = &(((CellStudyBWFFittingParameters*)data)->dySpectra);
	const std::vector<std::vector<double>>* const sf = &(((CellStudyBWFFittingParameters*)data)->survivingFraction);
	const std::vector<std::vector<double>>* const dose = &(((CellStudyBWFFittingParameters*)data)->dose);
	//const double* beta = &(((CellStudySurvivalParameters*)data)->beta);

	const BiologicalWeightingFunction* alphaFittingFunction = &(((CellStudyBWFFittingParameters*)data)->alphaFittingFunction);
	const BiologicalWeightingFunction* betaFittingFunction = &(((CellStudyBWFFittingParameters*)data)->betaFittingFunction);
	const int nAlphaParams = alphaFittingFunction->GetNumFittingParams(); //this is used to determine the offset in the parameter list, when computing beta
	const int nParams = (alphaFittingFunction->GetNumFittingParams()+betaFittingFunction->GetNumFittingParams()); 

	double params[nParams];

	//Pull the current fitting parameters
	for (int m = 0; m < nParams; ++m)
	{
		params[m] = gsl_vector_get(x,m);
	}

	//I know, this is confusing because I have 4 different iterators
	// i for the overall # of iterations, k for the histogram bins, j for the dose level, l for the LET
	int i = 0; //To keep track of which LET and Dose we are on
	int l = 0; //To keep track of just which LET we are on


	for(const std::pair<std::string,TH1D>& spectrumPair:*pDySpectra) //First we iterate over each of the lineal energy spectra
	{
		const TH1D& dySpectrum = std::get<1>(spectrumPair); //Pulling a reference so we don't remake this object
		double alphaPredicted = 0;
		double betaPredicted = 0;

		for(int k = 1; k <= dySpectrum.GetNbinsX(); ++k) //Iterate over every histogram bin to calculate alpha
		{
			//Get the spectrum value
			double width = dySpectrum.GetBinWidth(k);
			double center = dySpectrum.GetBinCenter(k);
			double value = dySpectrum.GetBinContent(k);

			//Accumulate the predicted value of alpha as we integrate the function
			double ryAlphaVal = alphaFittingFunction->GetValue(params,center);
			alphaPredicted += ryAlphaVal*value*width; //value*width is d(y)*dy, and everything else is r(y)

			double ryBetaVal = betaFittingFunction->GetValue(&(params[nAlphaParams]), center);
			betaPredicted += ryBetaVal*value*width;
		}

		for (int j = 0; j < (*dose)[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	
			double survivalPredicted = 0; 

			//calculate the exponent of survival predicted
			survivalPredicted = (alphaPredicted*(((*dose)[l])[j]))+((betaPredicted)*(((*dose)[l])[j])*(((*dose)[l])[j]));
			//take e^exponent
			survivalPredicted = std::exp(-survivalPredicted);

			//Add entry to cost function vector with value Y_Predicted - Y_GroundTruth
			gsl_vector_set(f,i,survivalPredicted-((*sf)[l])[j]); 
			++i;
		}

		++l; //iterate the LET
	}
	
	return GSL_SUCCESS;
}

double* BWF_Fitter_AlphaBeta::GeneralizedBWFFitting(CellStudyBWFFittingParameters& survivalParams, BiologicalWeightingFunction fitFuncAlpha, BiologicalWeightingFunction fitFuncBeta, double* initialGuess, bool weightedFit)
{
	survivalParams.alphaFittingFunction = fitFuncAlpha;
	survivalParams.betaFittingFunction = fitFuncBeta;

	//Setting up the fitting algorithm
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  	gsl_multifit_nlinear_workspace *w;
  	gsl_multifit_nlinear_fdf fdf;
 	gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

 	//Some tunable parameters
 	//fdf_params.factor_up = 1.5; fdf_params.factor_down = 1.5;

 	//Set tolerances
	const double xtol = 1e-8;
	const double gtol = 1e-8;
	const double ftol = 0.0;
    const size_t nDataPoints = survivalParams.GetNumDataPoints();
    int nFittingParams = (fitFuncAlpha.GetNumFittingParams()+fitFuncBeta.GetNumFittingParams());

  	//Set some initial values
  	int status, info;
  	double chisq, chisq0;
  	gsl_matrix *J;
  	gsl_matrix *covar = gsl_matrix_alloc (nFittingParams, nFittingParams);

	//Random number generator
	gsl_rng * r;
	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);

	//Set initial guess
	double* x_initial;
	if (initialGuess == nullptr)
	{
		x_initial = new double[nFittingParams];
		for (int i = 0; i < nFittingParams; ++i)
		{
			x_initial[i] = 1;
		}	
	}
	else
	{
		x_initial = initialGuess;
	}
	
	
	double weights[nDataPoints];
	gsl_vector_view x = gsl_vector_view_array (x_initial, nFittingParams);
    gsl_vector_view wts = gsl_vector_view_array(weights, nDataPoints);

	//Set weights if requested
	if (weightedFit)
	{
		int i = 0;
		for (const auto& experimentlist:survivalParams.survivingFractionUncertainty)
		{
			for (const auto& experimentuncertainty:experimentlist)
			{
				weights[i] = 1/(experimentuncertainty*experimentuncertainty);
				++i;
			}
		}
	} 
	else //otherwise weights are 1
	{
		for (int i = 0; i < nDataPoints; i++)
		{
			weights[i] = 1.0;
		}	
	}


	//Define function for minimization
	fdf.f = GeneralizedBWF;
	fdf.df = NULL;  //Not calculating the Jacobian explicitly
	fdf.fvv = NULL; //No geodesic acceleration
	fdf.n = nDataPoints; 
	fdf.p = nFittingParams; 
	fdf.params = &survivalParams;

	//allocate workspace with default parameters
  	w = gsl_multifit_nlinear_alloc (T, &fdf_params, nDataPoints, nFittingParams);

	//initialize solver with starting point and weights
  	gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

  	//compute initial cost function 
  	gsl_vector *f;
  	f = gsl_multifit_nlinear_residual(w);
  	gsl_blas_ddot(f, f, &chisq0);

  	//solve the system with a maximum of 100 iterations
  	status = gsl_multifit_nlinear_driver(400, xtol, gtol, ftol, callback, &nFittingParams, &info, w);

  	//Calculate the final cost
  	gsl_blas_ddot(f, f, &chisq);

  	//End of function stuff
  	fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(w),
          gsl_multifit_nlinear_trs_name(w));
	fprintf(stderr, "number of iterations: %zu\n",
	      gsl_multifit_nlinear_niter(w));
	fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
	fprintf(stderr, "reason for stopping: %s\n",
	      (info == 1) ? "small step size" : "small gradient");

	J = gsl_multifit_nlinear_jac(w);
 	gsl_multifit_nlinear_covar (J, 0.0, covar);
	#define FIT(i) gsl_vector_get(w->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	double* output = new double[nFittingParams];

	{
	    double dof = nDataPoints - nFittingParams;
	    double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

	    fprintf(stderr, "chisq/dof = %g\n", chisq / dof);
	    for (int i = 0; i < nFittingParams; ++i)
	    {
	    	fprintf (stdout, "C_%i     = %.5f +/- %.5f\n",i,FIT(i),c*ERR(i));
	    	output[i] = FIT(i);
	    }
	}
	
	fprintf (stderr, "status = %s\n", gsl_strerror (status));
	gsl_multifit_nlinear_free (w);
	gsl_rng_free (r);
	if (initialGuess == nullptr) {delete []initialGuess;} //Have to free the initial guesses if they weren't passed in

	return output;
}