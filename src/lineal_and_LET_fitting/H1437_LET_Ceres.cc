//std
#include <iostream>
#include <filesystem>
#include <vector>
#include <utility>
#include <string>
#include <cmath>

//ROOT
#include "TH1F.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TMath.h"

//This project
//Proton KESpectra
#include "ProtonSpectra.h" 
//Lineal energy spectra from KE Spectra
#include "LinealSpectra.h" 
//Functions for fitting and plotting
#include "SurvivalPlotting.h"
#include "Ceres_LET_Fitter.h"
#include "Ceres_LET_Fitter.h"
#include "AlphaBeta_Fitter.h"
#include "BWF_Fitting_Results.h"
#include "BWF_Plotter.h"
//Ceres Fitter and Google Logging
#include "ceres/ceres.h"
#include "glog/logging.h"
//This file
#include "H1437_LET_Ceres.h"

CellStudyBWFFittingParameters H1437FittingParams{};

void SetupH1437SurvivalParameters()
{
	//H1437 doses
	std::vector<std::vector<double>> doses;
	doses.push_back(std::vector<double>{0.4692676, 1.094958, 1.564225, 2.346338, 3.12845}); //5
	doses.push_back(std::vector<double>{0.3903299, 0.7806598, 1.366155, 1.951649, 2.927474, 3.903299}); //6
	doses.push_back(std::vector<double>{0.2615113,0.5230227,0.784534,1.046045,1.307557,1.569068,1.830579,2.092091,2.615113,3.92267,5.230227}); //11
	doses.push_back(std::vector<double>{0.292021,0.584042,0.876063,1.168084,1.460105,1.752126,2.044147,2.336168,2.92021,4.380315,5.84042}); //11
	doses.push_back(std::vector<double>{0.3109079,0.6218159,0.9327238,1.243632,1.55454,1.865448,2.176355,2.487263,3.109079,4.663619,6.218158}); //11
	doses.push_back(std::vector<double>{0.3699901,0.7399802,1.10997,1.47996,1.849951,2.219941,2.589931,2.959921,3.699901,5.549852,7.399802}); //11
	doses.push_back(std::vector<double>{0.4644248,0.9288495,1.393274,1.857699,2.322124,2.786549,3.250973,3.715398,4.644248,6.966372}); //10
	doses.push_back(std::vector<double>{0.6198787,1.239757,1.859636,2.479515,3.099394,3.719272,4.339151,4.95903,6.198787}); //9
	doses.push_back(std::vector<double>{0.8668616,1.733723,2.600585,3.467447,4.334308,5.20117,6.068032}); //7
	doses.push_back(std::vector<double>{0.687678,1.375356,2.063034,2.750712,3.43839,4.126068,4.813745}); //7
	doses.push_back(std::vector<double>{0.4281037,0.8562075,1.284311,1.712415,2.140519,2.568623,2.996726,3.42483}); //8
	doses.push_back(std::vector<double>{0.2619956,0.5239912,0.7859868,1.047982,1.309978,1.571974,1.833969,2.095965,2.619956}); //9

	//H1437 surviving fractions
	std::vector<std::vector<double>> survivingFractions;
	survivingFractions.push_back(std::vector<double>{0.98444715,0.89215520625,0.78960865625,0.74859000625,0.5906682625});
	survivingFractions.push_back(std::vector<double>{0.929071975,0.90856265625,0.787557725,0.67065459375,0.6296359625,0.41223721875});
	survivingFractions.push_back(std::vector<double>{1.0131601125,0.96803975,0.9844471625,0.8983080375,0.85728938125,0.886002475,0.81832165,0.73218256875,0.6768074,0.510681956875,0.346607388125});
	survivingFractions.push_back(std::vector<double>{1.06648441875,1.0787899375,0.92086821875,0.85113655,0.82037260625,0.77115024375,0.82242356875,0.6521962125,0.61117756875,0.37942229375,0.19894033125});
	survivingFractions.push_back(std::vector<double>{0.9598359875,1.0398223625,0.92291920625,0.84703475625,0.7752520875,0.7403862625,0.693214875,0.63989063125,0.5701589625,0.31789435,0.178431019375});
	survivingFractions.push_back(std::vector<double>{1.07879008125,0.94753036875,0.8880533625,0.8244744625,0.72808068125,0.7096223,0.61322851875,0.5927192,0.38352416875,0.231755231875,0.094342841875});
	survivingFractions.push_back(std::vector<double>{0.91881738125,0.89010426875,0.865493125,0.75679371875,0.57015895,0.58246455,0.52503846875,0.3609639,0.29943597,0.096393771875});
	survivingFractions.push_back(std::vector<double>{0.9905999,0.9126645375,0.838831025,0.600922925,0.52503845625,0.44710305625,0.3076396875,0.2297043,0.118954015});
	survivingFractions.push_back(std::vector<double>{0.92291916875,0.783455825,0.4696632875,0.3096906375,0.15792170875,0.149717985625,0.084088186875});
	survivingFractions.push_back(std::vector<double>{0.873696825,0.6932148375,0.42044093125,0.3076397,0.14561612375,0.1107502925,0.051273285625});
	survivingFractions.push_back(std::vector<double>{0.898308075,0.69936765,0.49837635,0.3773713625,0.2173987275,0.149717984375,0.108699363125,0.065629805});
	survivingFractions.push_back(std::vector<double>{0.902409875,0.869594975,0.5845154625,0.47581609375,0.34455648125,0.319945294375,0.2563664075,0.1210049475,0.05742608});


	//H1437 surviving fraction uncertainty (not yet filled out)
	// std::vector<std::vector<double>> survivingFractionUncertainty;
	// survivingFractionUncertainty.push_back(std::vector<double>{0.0973136077289898,0.103015986543154,0.103096420724027,0.107761376391554,0.0955386592610086,0.117579710048443,0.0695918896923402,0.0602326700663846,0.0489658382266086});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.109664062656694,0.114103383905432,0.106349029381286,0.100750909912505,0.104440896026642,0.115615500309136,0.0571078918837479,0.0596822436716957,0.0354580361405374});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.117538438828423,0.112939677206714,0.120312224123347,0.104189642227107,0.1110417317714,0.0699579930309851,0.058429879504437,0.0307981166308652,0.0143918780859903});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.116432603567846,0.0959351429808473,0.114867846060841,0.104554122926044,0.101609681173067,0.0806225179654492,0.0584744019992923,0.0206631372293671});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.107430799661687,0.112503375873942,0.118147015517618,0.101819298520923,0.0946038300253763,0.0919585883268045,0.06099555369411,0.0207280257711102});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.0977871206112248,0.100992226918609,0.115870692740975,0.0994505219983238,0.0929887914303216,0.0753478617253426,0.0313360366270534});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.120016044348618,0.112333967115251,0.0951047926938,0.0795554862324001,0.0768468409174345,0.0746030087767288,0.0246618969356351});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.117363566443651,0.105874866756647,0.0937622926169325,0.0809119467597451,0.0535464653407672,0.0418100763356018});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.0975556035780127,0.0760194642491641,0.055062426425359,0.0284367658829579,0.0190867076778482});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.10901048863062,0.0660380349795841,0.0342207978196262,0.0198752128828158});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.0756307934528414,0.0757123399282646,0.0522682804633686,0.0325580569599345,0.0162016671380295});
	// survivingFractionUncertainty.push_back(std::vector<double>{0.101115832072717,0.0867612613187053,0.0816322353302927,0.0437450743705149,0.0291949733842577,0.0224833348620257});


	//LETd for each experiment (from Fada's calculations)
	std::vector<double> LETd{0.914, 1.16, 1.61, 1.81, 1.93, 2.34, 3.01, 5.08, 10.8, 15.2, 17.7, 19};

	//These are actually the YD's I calculated
	// std::vector<double> LETd{1.92265, 2.15363, 2.53646,2.73784,2.89062,3.36408,4.17884,5.885,12.8591,18.2029,21.0286,22.7578};

	//Photon beta (from Fada's calculations)
	double beta = 0.097;

	//Manually fitted alpha-beta for each experiment (from Fada's calculations)
	std::vector<double> alphas {0.268, 0.226, 0.151, 0.150, 0.166, 0.137, 0.206, 0.117, 0.318, 0.446, 0.596, 0.883};
	std::vector<double> betas {0.097, 0.112, 0.134, 0.134, 0.134, 0.146, 0.125, 0.159, 0.154, 0.341, 0.662, 0.956};

	//Get the lineal energy spectra library
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");
	std::vector<std::pair<std::string, TH1D>>* dySpectraLibrary;
	TFile* dySpectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraCellStudy_1e3nm.root");
	dySpectrumFile->GetObject("Lineal_energy_library",dySpectraLibrary);
	dySpectrumFile->Close(); //Close file

	//Set the cell study parameters for the new fitting class too
	H1437FittingParams.dySpectra = *dySpectraLibrary;
	H1437FittingParams.dose = doses;
	H1437FittingParams.survivingFraction = survivingFractions;
	// H1437FittingParams.survivingFractionUncertainty = survivingFractionUncertainty;
	H1437FittingParams.beta = beta;	
	H1437FittingParams.LETd = LETd;
}

void Fixed_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	std::vector<double> quadraticGuesses;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.push_back(0);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Linear_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	std::vector<double> quadraticGuesses;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.push_back(0);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double>{0.0525633, -0.00443009, 0.0081919}); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	// LinearBWF.SetValues(std::vector<double> {0,0});
	// BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses{0.0768515, -0.000181136, 0.000284008, -0.00917733};
	LQE2BWF.SetValues(LQE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Quadratic_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	std::vector<double> quadraticGuesses;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.push_back(0);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses{0.103306, -0.145912, 0.00030309, -0.00372254};
	// LQEGuesses.push_back(0);
	// LQEGuesses.push_back(0.1);
	// LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	// LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double>{0.0805378, -0.00772915, 0.00224516}); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses{0.0768515, -0.000181136, 0.000284008, -0.00917733};
	LQE2BWF.SetValues(LQE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Cubic_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	std::vector<double> quadraticGuesses;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.push_back(0);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses{0.0884697, -0.00776365, 0.000113385, -0.000197189};
	LQE2BWF.SetValues(LQE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	CubicBWF.SetValues(std::vector<double> {0,0.000199399, 0.0002002, 0.0307028});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Fourth_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	// LinearBWF.SetValues(std::vector<double> {0.0129481, 0.0239448});
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FourthBWF.SetValues(std::vector<double> {0, 5.72865e-06, -0.000125374, -0.00368136, 0.050664});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	// std::vector<double> quadraticGuesses{0.00097301, -0.0267699, 0.132568};
	// QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	// FourthBWF.SetValues(std::vector<double> {0,6.54635e-06, -0.000487498, 0.00844529, 0.018007});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses{5.91184e-06, -0.0002469, 0.00344482, -0.0169155, 0.0205553, 0.104009};
	// fifthGuesses.push_back(0);
	// fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	// fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	// fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	// fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	// fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

	// 	//
	// 	// Mairani et al.
	// 	//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	// FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	// BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses{0.0984425, -0.011425, -3.79679e-06, 0.000457635};
	LQE2BWF.SetValues(LQE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// FourthBWF.SetValues(std::vector<double> {0,6.15596e-06, -0.000439335, 0.00710176, 0.0214961});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	// BWFFitter.SetPositiveConstrained(true, 0.01);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	// BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Fifth_Beta()
{
	SetupH1437SurvivalParameters();

	//
	// Creating BWFS
	//

		//
		// The Polynomials
		//

	//Create fixed "BWF"
	BiologicalWeightingFunction FixedBWF;
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);

	//Create linear BWF
	BiologicalWeightingFunction LinearBWF;
	LinearBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticBWF;
	QuadraticBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 3);

	//Create cubic BWF
	BiologicalWeightingFunction CubicBWF;
	CubicBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 4);

	//Fourth order poly BWF
	BiologicalWeightingFunction FourthBWF;
	FourthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 5);

	//Fifth order poly BWF
	BiologicalWeightingFunction FifthBWF;
	FifthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[5]*linealEnergy*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[4]*linealEnergy*linealEnergy*linealEnergy*linealEnergy)+(params[3]*linealEnergy*linealEnergy*linealEnergy)+(params[2]*linealEnergy*linealEnergy)+(params[1]*linealEnergy)+params[0]);}, 6);

		//
		// From Mairani et al. 2017
		//

	//Create Q BWF (just quadratic nothing else)
	BiologicalWeightingFunction QBWF;
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy+params[1]);}, 2);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy)+params[2]);}, 3);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy)+params[3]);}, 4);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy)+params[2]);}, 3);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy)+params[3]);}, 4);

		//
		// Other things I've tried
		//

	//Exponential growth BWF
	BiologicalWeightingFunction ExponentialGrowthBWF;
	ExponentialGrowthBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*std::exp(linealEnergy/params[1]));}, 2);

	//Exponential, starts at some value then decaying to zero BWF
	BiologicalWeightingFunction ExponentialToZeroBWF;
	ExponentialToZeroBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[1]*(1-std::exp(-params[0]/linealEnergy)) ;}, 2);

	//Simple exponential decay BWF
	BiologicalWeightingFunction ExponentialDecayBWF;
	ExponentialDecayBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return std::exp(-params[0]*linealEnergy) ;}, 1);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianBWF;
	GaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((1/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	}
	, 2);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	}
	, 3);

	//BiologicalWeightingFunction GaussianBWF;
	BiologicalWeightingFunction GaussianVariableAmplitudePlusOffsetBWF;
	GaussianVariableAmplitudePlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);

	BiologicalWeightingFunction GaussianPlusOffsetBWF;
	GaussianPlusOffsetBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0]/params[2])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))); 
	}
	, 3);

	//Three parameter sigmoid;
	BiologicalWeightingFunction ThreeParameterSigmoid;
	ThreeParameterSigmoid.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])/(1+std::exp(-params[1]*(linealEnergy-params[2])))); // a / ( 1+e(-b*(y-c)) )
	}
	, 3);

	BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy);
	}, 4);

	BiologicalWeightingFunction MorstinPlusOffsetBWF;
	MorstinPlusOffsetBWF.SetWeightingFunction([](double const* params, double linealEnergy) 
	{
		return (params[0]*((1-std::exp(-linealEnergy*params[1]-params[2]*linealEnergy*linealEnergy-params[3]*linealEnergy*linealEnergy*linealEnergy))/linealEnergy))+params[4];
	}, 5);

	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0];
	    double omega = params[1];
	    double alpha = params[2];
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	    return (2./omega) * smallphi * bigphi;
	}
	, 3);

	BiologicalWeightingFunction SkewGaussianVariableAmplitudeBWF;
	SkewGaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	    double xi = params[0]; //std::cout << "xi: " << xi << std::endl;
	    double omega = params[1]; ///std::cout << "omega: " << omega << std::endl;
 	    double alpha = params[2]; //std::cout << "alpha: " << alpha << std::endl;
	    double amplitude = params[3]; //std::cout << "amplitude: " << amplitude << std::endl;
	    double arg = (linealEnergy - xi) / omega;
	    double smallphi = std::exp(-((linealEnergy-xi)*(linealEnergy-xi))/(params[1]*params[1]*2)); 
	    // TMath::Gaus(arg, 0.0, 1.0, true);
	    double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));

	    return (2*amplitude) * smallphi * bigphi; //2* is so scale is the same as gaussian
	}
	, 4);

	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_LET_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	// LinearBWF.SetValues(std::vector<double> {0.0129481, 0.0239448});
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {0, 0, 5.72865e-06, -0.000125374, -0.00368136, 0.050664});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.SetCellStudyParameters(H1437FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	// std::vector<double> quadraticGuesses{0.00097301, -0.0267699, 0.132568};
	// QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 5.72865e-06, -0.000125374, -0.00368136, 0.050664});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QResults = BWFFitter.Fit();
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	std::vector<double> LEGuesses;
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	// FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	// BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	std::vector<double> QE2Guesses{0.0974585, -0.00844037, 5.60171e-05};
	QE2BWF.SetValues(QE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FifthBWF.SetValues(std::vector<double> {0,2.69047e-06, -2.81043e-05, -0.000293497, 0.0062276, 0.0178759});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses{0.116671, -8.06524e-05, 0.000625052, -0.0197886};
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// FifthBWF.SetValues(std::vector<double> {0,0,6.15596e-06, -0.000439335, 0.00710176, 0.0214961});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 105);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_LET_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

	// std::cout << "Skew Gaussian" << std::endl;
	// SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {6, -5, 35, 100}); //SkewGaussianVariableAmplitudeBWF.SetValues(std::vector<double> {3, 0, 40, 80});  initial guess that does someething
	// BWFFitter.SetAlphaWeightingFunction(SkewGaussianVariableAmplitudeBWF);
	// // FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 1, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(3,1000); 	//Prevent runaway amplitude
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto SkewGaussianVariableAmplitudeResults = BWFFitter.Fit();
	// SkewGaussianVariableAmplitudeResults.PrintBasic();

	//
	// Plotting
	//

		//
		// Plotting the Polynomial BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//Plotting Beta

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c2","c2");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// lineStyle.SetLineColor(kGreen+2);
	// // BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// // lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, cubicResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, linearResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, quadraticResults.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialBeta_fifthAlpha_BetaBWFs.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting the Mairani BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.alphaFunc, QResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.alphaFunc, LEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.alphaFunc, QEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.alphaFunc, LQEResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.alphaFunc, LE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.alphaFunc, QE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.alphaFunc, LQE2Results.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/AlphaMairani_BetaFifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Beta
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c2 = new TCanvas("c","c");
	// c2->SetCanvasSize(9000, 5000);
	// c2->SetFillStyle(4000);
	// c2->SetFrameFillStyle(4000);
	// c2->Divide(4,3,0.000000005,0.001);

	// BWFFunctionPlotter(c, legend, lineStyle, "", QResults.betaFunc, QResults.betaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LEResults.betaFunc, LEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QEResults.betaFunc, QEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQEResults.betaFunc, LQEResults.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kYellow+1);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LE2Results.betaFunc, LE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kAzure+10);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QE2Results.betaFunc, QE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kBlue-8);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LQE2Results.betaFunc, LQE2Results.betaFunc.GetFittingParams(), "L", 0., 100.);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BetaFifth_AlphaMairani.jpg";
	// c2->SaveAs((TString)outputName); 

		//
		// Plotting other BWFs
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinResults.alphaFunc, MorstinResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Morstin_alpha_to_250.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.betaFunc, GaussianVarAmplitudeResults.betaFunc.GetFittingParams(), "AL", 0.01, 110.);	
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_beta.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", MorstinOffsetResults.alphaFunc, MorstinOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/mostinoffset_positive_constrained.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudePlusOffsetResults.alphaFunc, GaussianVarAmplitudePlusOffsetResults.alphaFunc.GetFittingParams(), "AL", 0.01, 250.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_offset_fifth_bwf.jpg";
	// c->SaveAs((TString)outputName); 

	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianVariableAmplitudeResults.alphaFunc, SkewGaussianVariableAmplitudeResults.alphaFunc.GetFittingParams(), "AL", 0.01, 100.);	
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVarAmplitudeResults.alphaFunc, GaussianVarAmplitudeResults.alphaFunc.GetFittingParams(), "L", 0.01, 100.);	
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/skew_gaussian_fifth_BWF_0to100.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the Polynomial Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H1437FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H1437FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H1437FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H1437FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H1437FittingParams, fifthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	// c->SaveAs((TString)outputName); 


		//
		// Plotting the Mairani Residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H1437FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H1437FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H1437FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H1437FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H1437FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H1437FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H1437FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting the other residuals
		//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	// TAttLine lineStyle{};
	// lineStyle.SetLineColor(kGreen+2);

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H1437FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H1437FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H1437FittingParams);
	// double* AlphaBeta = fitter.Fit(nullptr,false);

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);
	// TAttMarker markerAtts;
	// markerAtts.SetMarkerColor(kBlack);
	// markerAtts.SetMarkerSize(8);
	// markerAtts.SetMarkerStyle(8);

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H1437FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H1437FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H1437FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void H1437_LET_Ceres()
{
	Fifth_Beta();
}