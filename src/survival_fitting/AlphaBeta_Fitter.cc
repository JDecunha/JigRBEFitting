#include "AlphaBeta_Fitter.h"

double* AlphaBeta_Fitter::Fit(double* initialGuess, bool weightedFit)
{
	if(!_paramsSet) {std::cout << "Attempt to use BWF_Fitter_AlphaBeta without setting survival params. Failure." << std::endl; return nullptr; }

	return LQModelFitting(_survivalParams, initialGuess, weightedFit);
}

int AlphaBeta_Fitter::LQModel(const gsl_vector* x, void* data, gsl_vector* f)
{
	//Load in the surviving fraction, dose, and LET
	const std::vector<std::vector<double>>* const sf = &(((CellStudyBWFFittingParameters*)data)->survivingFraction);
	const std::vector<std::vector<double>>* const dose = &(((CellStudyBWFFittingParameters*)data)->dose);
	const std::vector<double>* const LETds = &(((CellStudyBWFFittingParameters*)data)->LETd);
	const int nParams = ((*LETds).size()*2);

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


	for(double LETd:*LETds) //First we iterate over each LET
	{

		for (int j = 0; j < (*dose)[l].size(); ++j) //Iterate over every different dose level and determine surviving faction
		{	
			double survivalPredicted = 0; 
			double alpha = params[2*l]; //starts at 0, increments by 2
			double beta = params[2*l+1]; //starts at 1, increments by 1

			//calculate the exponent of survival predicted
			survivalPredicted = (alpha*(((*dose)[l])[j]))+((beta)*(((*dose)[l])[j])*(((*dose)[l])[j]));
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

double* AlphaBeta_Fitter::LQModelFitting(CellStudyBWFFittingParameters& survivalParams, double* initialGuess, bool weightedFit)
{

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
    std::vector<double> LETds = survivalParams.LETd;
    int nFittingParams = ((LETds).size()*2);

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
	fdf.f = LQModel;
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