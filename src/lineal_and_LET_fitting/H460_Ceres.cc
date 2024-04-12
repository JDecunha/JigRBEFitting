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
#include "Ceres_Survival_Fitter.h"
#include "Ceres_BWF_Fitter.h"
#include "Ceres_LET_Fitter.h"
#include "AlphaBeta_Fitter.h"
#include "BWF_Fitting_Results.h"
#include "BWF_Plotter.h"
//Ceres Fitter and Google Logging
#include "ceres/ceres.h"
#include "glog/logging.h"
//This file
#include "H460_Ceres.h"

CellStudyBWFFittingParameters H460FittingParams{};

void SetupH460SurvivalParameters()
{
	//H460 doses
	std::vector<std::vector<double>> doses;
	doses.push_back(std::vector<double>{0.1564225,0.3128451,0.4692676,0.6256901,0.7821126,0.9385352,1.564225,2.346338,3.12845}); //9
	doses.push_back(std::vector<double>{0.1951649,0.3903299,0.5854948,0.7806598,0.9758247,1.17099,1.951649,2.927474,3.903299}); //9
	doses.push_back(std::vector<double>{0.2615113,0.5230227,0.784534,1.046045,1.307557,1.569068,2.615113,3.92267,5.230227}); //9
	doses.push_back(std::vector<double>{0.292021,0.584042,0.876063,1.168084,1.460105,1.752126,2.92021,4.380315}); //8
	doses.push_back(std::vector<double>{0.3109079,0.6218159,0.9327238,1.243632,1.55454,1.865448,3.109079,4.663619}); //8
	doses.push_back(std::vector<double>{0.3699901,0.7399802,1.10997,1.47996,1.849951,2.219941,3.699901}); //7
	doses.push_back(std::vector<double>{0.4644248,0.9288495,1.393274,1.857699,2.322124,2.786549,4.644248}); //7
	doses.push_back(std::vector<double>{0.6198787,1.239757,1.859636,2.479515,3.099394,3.719272}); //6
	doses.push_back(std::vector<double>{0.8668616,1.733723,2.600585,3.467447,4.334308}); //5
	doses.push_back(std::vector<double>{0.687678,1.375356,2.063034,2.750712}); //4
	doses.push_back(std::vector<double>{0.4281037,0.8562075,1.284311,1.712415,2.140519}); //5
	doses.push_back(std::vector<double>{0.2619956,0.5239912,0.7859868,1.047982,1.309978,1.571974}); //6
	
	//H460 surviving fractions
	std::vector<std::vector<double>> survivingFractions;
	survivingFractions.push_back(std::vector<double>{0.93307278125,0.8790838890625,0.84050740625,0.8517369109375,0.7641728796875,0.721009665625,0.51075015625,0.324695915625,0.160247745});
	survivingFractions.push_back(std::vector<double>{0.9100104015625,0.902505121875,0.8416153125,0.7996460109375,0.707484875,0.680355915625,0.415398284375,0.201869659375,0.0734918634375});
	survivingFractions.push_back(std::vector<double>{0.95365396875,0.914083753125,0.8141567421875,0.7670349109375,0.6474413015625,0.522880384375,0.263363778125,0.0662321278125,0.0166767615});
	survivingFractions.push_back(std::vector<double>{0.93415940625,0.884765775,0.8070793625,0.7203413078125,0.5987456484375,0.456758728125,0.207980753125,0.0432527975});
	survivingFractions.push_back(std::vector<double>{0.9192694609375,0.8731507953125,0.781966409375,0.6774832546875,0.5528024296875,0.407313846875,0.1659575775,0.028806982125});
	survivingFractions.push_back(std::vector<double>{0.9364116578125,0.831354428125,0.7396655140625,0.6250709890625,0.45639296875,0.2882206,0.1010323921875});
	survivingFractions.push_back(std::vector<double>{0.8952672015625,0.749569390625,0.609157578125,0.4515768609375,0.29223910625,0.187674129375,0.03900813625});
	survivingFractions.push_back(std::vector<double>{0.856657646875,0.6869433875,0.45702575625,0.2987832265625,0.1351156596875,0.078734106875});
	survivingFractions.push_back(std::vector<double>{0.7073390125,0.3618327734375,0.134352829375,0.0540511246875,0.0193577265625});
	survivingFractions.push_back(std::vector<double>{0.646620878125,0.270156696875,0.0917367021875,0.02544827584375});
	survivingFractions.push_back(std::vector<double>{0.6958057296875,0.3780136796875,0.13946794375,0.05421250125,0.01603564625});
	survivingFractions.push_back(std::vector<double>{0.7483156359375,0.487375509375,0.2871622453125,0.1245501896875,0.055777303125,0.0375871184375});

	//H460 surviving fraction uncertainty
	std::vector<std::vector<double>> survivingFractionUncertainty;
	survivingFractionUncertainty.push_back(std::vector<double>{0.0973136077289898,0.103015986543154,0.103096420724027,0.107761376391554,0.0955386592610086,0.117579710048443,0.0695918896923402,0.0602326700663846,0.0489658382266086});
	survivingFractionUncertainty.push_back(std::vector<double>{0.109664062656694,0.114103383905432,0.106349029381286,0.100750909912505,0.104440896026642,0.115615500309136,0.0571078918837479,0.0596822436716957,0.0354580361405374});
	survivingFractionUncertainty.push_back(std::vector<double>{0.117538438828423,0.112939677206714,0.120312224123347,0.104189642227107,0.1110417317714,0.0699579930309851,0.058429879504437,0.0307981166308652,0.0143918780859903});
	survivingFractionUncertainty.push_back(std::vector<double>{0.116432603567846,0.0959351429808473,0.114867846060841,0.104554122926044,0.101609681173067,0.0806225179654492,0.0584744019992923,0.0206631372293671});
	survivingFractionUncertainty.push_back(std::vector<double>{0.107430799661687,0.112503375873942,0.118147015517618,0.101819298520923,0.0946038300253763,0.0919585883268045,0.06099555369411,0.0207280257711102});
	survivingFractionUncertainty.push_back(std::vector<double>{0.0977871206112248,0.100992226918609,0.115870692740975,0.0994505219983238,0.0929887914303216,0.0753478617253426,0.0313360366270534});
	survivingFractionUncertainty.push_back(std::vector<double>{0.120016044348618,0.112333967115251,0.0951047926938,0.0795554862324001,0.0768468409174345,0.0746030087767288,0.0246618969356351});
	survivingFractionUncertainty.push_back(std::vector<double>{0.117363566443651,0.105874866756647,0.0937622926169325,0.0809119467597451,0.0535464653407672,0.0418100763356018});
	survivingFractionUncertainty.push_back(std::vector<double>{0.0975556035780127,0.0760194642491641,0.055062426425359,0.0284367658829579,0.0190867076778482});
	survivingFractionUncertainty.push_back(std::vector<double>{0.10901048863062,0.0660380349795841,0.0342207978196262,0.0198752128828158});
	survivingFractionUncertainty.push_back(std::vector<double>{0.0756307934528414,0.0757123399282646,0.0522682804633686,0.0325580569599345,0.0162016671380295});
	survivingFractionUncertainty.push_back(std::vector<double>{0.101115832072717,0.0867612613187053,0.0816322353302927,0.0437450743705149,0.0291949733842577,0.0224833348620257});


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
	TFile* dySpectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraApril2023CellStudy_1e3nm.root");
	dySpectrumFile->GetObject("Lineal_energy_library",dySpectraLibrary);
	dySpectrumFile->Close(); //Close file

	//Set the cell study parameters for the new fitting class too
	H460FittingParams.dySpectra = *dySpectraLibrary;
	H460FittingParams.dose = doses;
	H460FittingParams.survivingFraction = survivingFractions;
	H460FittingParams.survivingFractionUncertainty = survivingFractionUncertainty;
	H460FittingParams.beta = beta;	
	H460FittingParams.LETd = LETd;
}

void Fixed_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	// LinearBWF.SetValues(std::vector<double> {0.0399,0.0457});
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.1355}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

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
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

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
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

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
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FixedBWF.SetValues(std::vector<double> {0});
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Linear_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

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
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

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
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// LinearBWF.SetValues(std::vector<double> {0,0});
	// BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

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
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	LinearBWF.SetValues(std::vector<double> {0,0});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Quadratic_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.01; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

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
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize(); 
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	// std::vector<double> LQE2Guesses;
	// LQE2Guesses.push_back(0.264955);
	// LQE2Guesses.push_back(-5.87335e-05);
	// LQE2Guesses.push_back(0.00166188);
	// LQE2Guesses.push_back(-0.0441215);
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// QuadraticBWF.SetValues(std::vector<double> {0,0.00455537, 0.124748});
	// BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

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
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,0,0});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Cubic_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.5; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

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
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

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
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	std::vector<double> QE2Guesses;
	QE2Guesses.push_back(.148733);
	QE2Guesses.push_back(-0.000220853);
	QE2Guesses.push_back(0.00014253);
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	CubicBWF.SetValues(std::vector<double> {0,0.000718521, -0.0217446, 0.19776});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	std::vector<double> LQE2Guesses;
	LQE2Guesses.push_back(0.232922);
	LQE2Guesses.push_back(-0.000114764);
	LQE2Guesses.push_back(0.00103702);
	LQE2Guesses.push_back(-0.0311275);
	LQE2BWF.SetValues(LQE2Guesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	CubicBWF.SetValues(std::vector<double> {0, 0.000361983, -0.00654466, 0.154374});
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

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
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);   
}

void Fourth_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	LinearBWF.SetValues(std::vector<double> {0.0250563, 0.0880142}); //values from cubic beta fitting
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FourthBWF.SetValues(std::vector<double> {0,3.94134e-06, 0.00110483, -0.0339738, 0.220536}); //values from cubic beta
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.1; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

		//
		// Mairani et al.
		//

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

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
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {0.150172, -0.000226012, 0.000121788});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FourthBWF.SetValues(std::vector<double> {0,8.29326e-06, 0.000106663, -0.0147373, 0.183807});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {0.228648, -0.000117843, 0.000980533, -0.0299326});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FourthBWF.SetValues(std::vector<double> {0,1.37078e-05, -0.00054164, -0.000957938, 0.150888});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto MorstinResults = BWFFitter.Fit();
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FourthBWF.SetValues(std::vector<double> {0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Fifth_Beta()
{
	SetupH460SurvivalParameters();

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

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear" << std::endl;
	LinearBWF.SetValues(std::vector<double> {0.0250563, 0.0880142}); //values from cubic beta fitting
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,3.94134e-06, 0.00110483, -0.0339738, 0.220536}); //values from cubic beta
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.1; //0.01;
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

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
	quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

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
	cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

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
	fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

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
	fifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthResults.betaFunc);

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
	QResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QResults.betaFunc);

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
	LEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEResults.betaFunc);

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
	LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2Results.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {0.150172, -0.000226012, 0.000121788});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FifthBWF.SetValues(std::vector<double> {0,0,8.29326e-06, 0.000106663, -0.0147373, 0.183807});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2Results.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {0.228648, -0.000117843, 0.000980533, -0.0299326});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

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
	MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FifthBWF.SetValues(std::vector<double> {0,0,0,0,0,0});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetResults.betaFunc);

	// std::cout << std::endl << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// BWFFitter.Initialize();
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianGaussianVarAmplitudeResults.betaFunc);


	// std::cout << std::endl << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudeResults.betaFunc);

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

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	legend->SetTextSize(0.05);

	TAttLine lineStyle{};
	lineStyle.SetLineColor(kGreen+2);

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Linear", H460FittingParams, linearResults);
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fifth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fifth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fifth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fifth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fifth.jpg";
	c->SaveAs((TString)outputName); 


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

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fifth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
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

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinResults);
	// // std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstin_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// // AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Morstin", H460FittingParams, MorstinOffsetResults);
	// // outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_morstinoffset_fifth.jpg";
	// // c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Gaussian", H460FittingParams, GaussianVarAmplitudeResults);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_gaussian_fifth.jpg";
	// c->SaveAs((TString)outputName); 

		//
		// Plotting Alpha and Beta Compared to direct fitting
		//

	// AlphaBeta_Fitter fitter{};
	// fitter.SetCellStudyParameters(H460FittingParams);
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

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, false);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_alpha.jpg";
	// c->SaveAs((TString)outputName); 

	// c->Clear();
	// legend->Clear();

	// markerAtts.SetMarkerColor(kBlack);
	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// markerAtts.SetMarkerColor(kRed+2);
	// PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, GaussianVarAmplitudeResults, true);
	// // markerAtts.SetMarkerColor(kRed+2);
	// // PlotAlphaBeta(c, legend, "", markerAtts, "P", H460FittingParams, AlphaBeta, true);
	// // markerAtts.SetMarkerColor(kGreen+2);
	// // PlotAlphaBetaFromBWF(c, legend, "", markerAtts, "P", H460FittingParams, MorstinResults, true);

	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/gaussian_fifth_beta.jpg";
	// c->SaveAs((TString)outputName);  
}

void Survival_Fitting()
{
	SetupH460SurvivalParameters();

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

	AlphaBeta_Fitter fitter{};
	fitter.SetCellStudyParameters(H460FittingParams);
	double* AlphaBeta = fitter.Fit(nullptr,false);

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->SetBottomMargin(0.13);
	// c->Divide(4,3,0.000000005,0.001);
	// auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	auto legend = new TLegend(0.14,0.64,0.14+0.31,0.72+0.23);
	legend->SetTextSize(0.04);

	TAttMarker markerAtts;

	markerAtts.SetMarkerSize(10);
	
	TAttLine lineStyle{};
	lineStyle.SetLineWidth(5);
	lineStyle.SetLineStyle(1);

	Ceres_Survival_Fitter CeresSurvival{};
	CeresSurvival.SetCellStudyParameters(H460FittingParams);
	CeresSurvival.Initialize();
	// auto results = CeresSurvival.Fit();

	// SurvivalDataMultigraph(c, legend, H460FittingParams);
	// lineStyle.SetLineColor(kRed+2);
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// MultigraphSurvivalFitPlotter(c, legend, lineStyle, "DirectFit", H460FittingParams, AlphaBeta);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/DirectFitH460.jpg";
	// c->SaveAs((TString)outputName); 

	// results.PrintSummary();
	

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/alpha_H460.jpg";
	// c->SaveAs((TString)outputName); 

	// PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/beta_H460.jpg";
	// c->SaveAs((TString)outputName); 

	double* cesiumAlphaBeta = new double[2];
	cesiumAlphaBeta[0] = 0.290; cesiumAlphaBeta[1] = 0.083;



	Ceres_BWF_Fitter BWFFitter{};
	BiologicalWeightingFunction GaussianPlusOffsetBWF2;
	BWF_Fitting_Results GaussianDyBestResults;
	GaussianPlusOffsetBWF2.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);
	GaussianPlusOffsetBWF2.SetValues(std::vector<double>  {  0.230491501205610, 5.51962016177844,114.092470171093,95.2590461228117 });
	FixedBWF.SetValues(std::vector<double> {0.101961735572583});
	GaussianDyBestResults.alphaFunc = GaussianPlusOffsetBWF2;
	GaussianDyBestResults.betaFunc = FixedBWF;

	BWF_Fitting_Results MorstinLinearDyBestResults;
	MorstinBWF.SetValues(std::vector<double> {1.33638439092475e-06, -2.98606802538980e-05, 0.000137006580329845, 3010.90150042332});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.111377434741229});
	MorstinLinearDyBestResults.alphaFunc = MorstinBWF;
	MorstinLinearDyBestResults.betaFunc = FixedBWF;

	BWF_Fitting_Results MorstinLinearDyBestWeightedResults;
	MorstinBWF.SetValues(std::vector<double> {7.87804976156578e-07, -2.13716006297374e-05, 0.000114200443445769, 5424.35985122075});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.0893680456034658});
	MorstinLinearDyBestWeightedResults.alphaFunc = MorstinBWF;
	MorstinLinearDyBestWeightedResults.betaFunc = FixedBWF;

	BWF_Fitting_Results CubicLinearDyBestRegularResults;
	// CubicBWF.SetValues(std::vector<double> {0.00838701932847768,0.185289234216632,-0.0250980413297348,0.000372587152280048});
	CubicBWF.SetValues(std::vector<double> {0.000337133584729527, -0.0222918654682546, 0.157099764860757, 0.0431408242731182});
	//FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.112398569166481});
	CubicLinearDyBestRegularResults.alphaFunc = CubicBWF;
	CubicLinearDyBestRegularResults.betaFunc = FixedBWF;

	BWF_Fitting_Results CubicLinearDyBestWeightedResults;
	// CubicBWF.SetValues(std::vector<double> {0.00838701932847768,0.185289234216632,-0.0250980413297348,0.000372587152280048});
	CubicBWF.SetValues(std::vector<double> {0.000372587152280048, -0.0250980413297348, 0.185289234216632, 0.00838701932847768});
	//FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.115994363406138});
	CubicLinearDyBestWeightedResults.alphaFunc = CubicBWF;
	CubicLinearDyBestWeightedResults.betaFunc = FixedBWF;

	Ceres_LET_Fitter BWFFitter2;
	QEBWF.SetValues(std::vector<double>  {  0.209024, -0.105069, 0.000636537 });
	FixedBWF.SetValues(std::vector<double> {0.114731});
	BWF_Fitting_Results BestLETH460;
	BestLETH460.alphaFunc = QEBWF;
	BestLETH460.betaFunc = FixedBWF;

	// markerAtts.SetMarkerColor(kGreen+2);
	lineStyle.SetLineWidth(14);
	gStyle->SetLineStyleString(11,"400 100");
	lineStyle.SetLineStyle(11);
	lineStyle.SetLineColor(kRed-6);
	PlotRBE10SFMcNamara(c, legend, "McNamara", markerAtts, lineStyle, "AL", H460FittingParams,cesiumAlphaBeta);
	lineStyle.SetLineColor(kAzure-8);
	PlotRBE10SFChenAndAhmad(c, legend, "Chen and Ahmad", markerAtts, lineStyle, "L", H460FittingParams,cesiumAlphaBeta);
	lineStyle.SetLineColor(kGreen-8);
	PlotRBE10SFWedenberg(c, legend, "Wedenberg", markerAtts, lineStyle, "L", H460FittingParams,cesiumAlphaBeta);

	markerAtts.SetMarkerColor(kAzure+2);
	markerAtts.SetMarkerStyle(323);
	markerAtts.SetMarkerSize(12);
	PlotRBE10SF(c, legend, "d(y) - Morstin, weighted", markerAtts, "P", H460FittingParams,cesiumAlphaBeta,MorstinLinearDyBestWeightedResults);

	markerAtts.SetMarkerColor(kTeal+4);
	markerAtts.SetMarkerStyle(34);
	markerAtts.SetMarkerSize(10);
	PlotRBE10SFLET(c, legend, "LET_{d}", markerAtts, "P", H460FittingParams,cesiumAlphaBeta,BestLETH460);
	// PlotRBE10SF(c, legend, "d(y) - Morstin", markerAtts, "P", H460FittingParams,cesiumAlphaBeta,MorstinLinearDyBestResults);
	markerAtts.SetMarkerColor(kBlack);markerAtts.SetMarkerStyle(23);
	PlotRBE10SF(c, legend, "Guan 2015", markerAtts, lineStyle, "P", H460FittingParams,cesiumAlphaBeta,AlphaBeta);
	

	

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/RBE_10_SF_H460_MorstinWeightedVsNon.jpg";
	// c->SaveAs((TString)outputName); 

	// void PlotRBE10SF(TCanvas* c, TLegend* legend, std::string const& legendName, TAttMarker const& markerAttributes, std::string options, CellStudyBWFFittingParameters survivalParams, double const* alphaBetasCesium, double const* alphaBetasProton)''
}

void Skew_Gaussian_Testing()
{
	BiologicalWeightingFunction SkewGaussianBWF;
	SkewGaussianBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
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

	BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 

	}
	, 3);

	// BiologicalWeightingFunction GaussianVariableAmplitudeBWF;
	// GaussianVariableAmplitudeBWF.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	// { 
	// 	// return ((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2))); 
	// 	return ((params[0])*std::exp(-(1./2.)*((linealEnergy-params[1])/params[2])*((linealEnergy-params[1])/params[2])));
	// }
	// , 3);


	//Manually set BWF values
	GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {88.1875, 372.485, 1000});
	SkewGaussianBWF.SetValues(std::vector<double> {8, -15, 35, 130}); //Param[0] 

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	legend->SetTextSize(0.05);

	TAttLine lineStyle{};
	lineStyle.SetLineWidth(5);
	lineStyle.SetLineStyle(1);


	lineStyle.SetLineColor(kGreen+2);
	BWFFunctionPlotter(c, legend, lineStyle, "", SkewGaussianBWF, SkewGaussianBWF.GetFittingParams(), "AL", 0.01, 250.);	
	lineStyle.SetLineColor(kPink+2);
	lineStyle.SetLineWidth(12);
	BWFFunctionPlotter(c, legend, lineStyle, "", GaussianVariableAmplitudeBWF, GaussianVariableAmplitudeBWF.GetFittingParams(), "L", 0.01, 250.);	
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/testing_gaussian_fits.jpg";
	c->SaveAs((TString)outputName); 
}

void Best_Model_Plotting()
{

	SetupH460SurvivalParameters();

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

	// Ceres_BWF_Fitter BWFFitter{};
	// std::cout << std::endl << "Linear" << std::endl;
	// BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// // FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// // FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	// BWFFitter.SetCellStudyParameters(H460FittingParams);
	// double penaltyWeight = 0.001; //0.01;
	// BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	// BWFFitter.Initialize();
	// auto linearResults = BWFFitter.Fit();
	// linearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(linearResults.betaFunc);

	// std::cout << std::endl << "Quadratic" << std::endl;
	// std::vector<double> quadraticGuesses;
	// //Take the results from the last fit and use them to inform this fit
	// quadraticGuesses.push_back(0);
	// quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	// quadraticGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[0]);
	// QuadraticBWF.SetValues(quadraticGuesses);
	// BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	// BWFFitter.Initialize();
	// auto quadraticResults = BWFFitter.Fit();
	// quadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResults.betaFunc);

	// std::cout << std::endl << "Cubic" << std::endl;
	// std::vector<double> cubicGuesses;
	// cubicGuesses.push_back(0);
	// cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	// cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[1]);
	// cubicGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[0]);
	// CubicBWF.SetValues(cubicGuesses);
	// BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	// BWFFitter.Initialize();
	// auto cubicResults = BWFFitter.Fit();
	// cubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(cubicResults.betaFunc);

	// std::cout << std::endl << "Fourth order polynomial" << std::endl;
	// std::vector<double> fourthGuesses;
	// fourthGuesses.push_back(0);
	// fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[3]);
	// fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[2]);
	// fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[1]);
	// fourthGuesses.push_back(cubicResults.alphaFunc.GetFittingParams()[0]);
	// FourthBWF.SetValues(fourthGuesses);
	// BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	// BWFFitter.Initialize();
	// auto fourthResults = BWFFitter.Fit();
	// fourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(fourthResults.betaFunc);

	// 	//
	// 	// Mairani et al.
	// 	//

	// std::cout << std::endl << "LQE" << std::endl;
	// std::vector<double> LQEGuesses;
	// LQEGuesses.push_back(0);
	// LQEGuesses.push_back(0.1);
	// LQEGuesses.push_back(quadraticResults.alphaFunc.GetFittingParams()[2]);
	// LQEGuesses.push_back(linearResults.alphaFunc.GetFittingParams()[1]);
	// LQEBWF.SetValues(LQEGuesses);
	// BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	// BWFFitter.Initialize(); 
	// auto LQEResults = BWFFitter.Fit();
	// LQEResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(LQEResults.betaFunc);

	// std::cout << std::endl << "LQE2" << std::endl;
	// LQE2BWF.SetValues(LQEGuesses);
	// BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	// BWFFitter.Initialize();
	// auto LQE2Results = BWFFitter.Fit();
	// LQE2Results.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2Results.betaFunc);

	//
	// Other
	//


	std::cout << std::endl << "Morstin" << std::endl;
	// MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	MorstinBWF.SetValues(std::vector<double> {1.33638439092475e-06, -2.98606802538980e-05, 0.000137006580329845, 3010.90150042332}); 
	// BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	CubicBWF.SetValues(std::vector<double> {0,0,0,0});
	// BWFFitter.SetBetaWeightingFunction(CubicBWF);
	// BWFFitter.Initialize();
	// auto MorstinResults = BWFFitter.Fit();
	// MorstinResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinResults.betaFunc);

	//1st order beta fitting

	// BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// // FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	// // FifthBWF.SetValues(std::vector<double> {-3.62923e-21, -1.02135e-05, 0.00129381, -0.0357744, 0.220213, -0.107979});
	// BWFFitter.SetBetaWeightingFunction(LinearBWF);
	// BWFFitter.Initialize();
	// auto linearResultsfirst = BWFFitter.Fit();
	// linearResultsfirst.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(linearResultsfirst.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(linearResultsfirst.betaFunc);

	// std::cout << std::endl << "Quadratic" << std::endl;
	// std::vector<double> quadraticGuessesfirst;
	// //Take the results from the last fit and use them to inform this fit
	// quadraticGuessesfirst.push_back(0);
	// quadraticGuessesfirst.push_back(linearResultsfirst.alphaFunc.GetFittingParams()[1]);
	// quadraticGuessesfirst.push_back(linearResultsfirst.alphaFunc.GetFittingParams()[0]);
	// QuadraticBWF.SetValues(quadraticGuessesfirst);
	// BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	// BWFFitter.Initialize();
	// auto quadraticResultsfirst = BWFFitter.Fit();
	// quadraticResultsfirst.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResultsfirst.alphaFunc);
	// Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticResultsfirst.betaFunc);

	BiologicalWeightingFunction GaussianPlusOffsetBWF2;
	GaussianPlusOffsetBWF2.SetWeightingFunction( [] (double const* params, double linealEnergy) 
	{ 
		return (((params[0])*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[2]*params[2]*2)))+params[3]); 
	}
	, 4);
	GaussianPlusOffsetBWF2.SetValues(std::vector<double> {  0.230491501205610, 5.51962016177844,114.092470171093,95.2590461228117 });

	//
	// Plotting
	//

		//
		// Plotting the Alpha BWFs
		//

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->SetBottomMargin(0.13);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// auto legend = new TLegend(0.14,0.66,0.14+0.33,0.72+0.22);
	legend->SetTextSize(0.06);

	TAttLine lineStyle{};
	lineStyle.SetLineWidth(5);
	lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "Quadratic", QuadraticBWF, quadraticResultsfirst.alphaFunc.GetFittingParams(), "AL", 0.001, 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "3rd", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "4th", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "LQE", LQEBWF, LQEResults.alphaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "LQE^{2}", LQE2BWF, LQE2Results.alphaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kYellow+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", GaussianPlusOffsetBWF2, GaussianPlusOffsetBWF2.GetFittingParams(), "AL", 0.001, 250.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Gaussian_test.jpg";
	// c->SaveAs((TString)outputName); 

	// TAttLine lineStyle{};
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);
	// lineStyle.SetLineColor(kYellow+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "Morstin", CubicBWF, MorstinResults.betaFunc.GetFittingParams(), "AL", 0.001, 100.);
	// lineStyle.SetLineColor(kGreen+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "Quadratic", LinearBWF, quadraticResultsfirst.betaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kPink+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "3rd", FixedBWF, cubicResults.betaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "4th", FixedBWF, fourthResults.betaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "LQE", FixedBWF, LQEResults.betaFunc.GetFittingParams(), "L", 0.001, 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "LQE^{2}", FixedBWF, LQE2Results.betaFunc.GetFittingParams(), "L", 0.001, 100.);
	

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/H460_Candidate_BetaBWFs.jpg";
	// c->SaveAs((TString)outputName); 


	//
	// Residuals
	//

	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(9000, 5000);
	// c->SetFillStyle(4000);
	// c->SetFrameFillStyle(4000);
	// c->Divide(4,3,0.000000005,0.001);
	// // auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// // auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	// legend->SetTextSize(0.05);

	

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "", H460FittingParams, quadraticResultsfirst);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, MorstinResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Morstin.jpg";
	// c->SaveAs((TString)outputName); 

	//
	//Survival data directly
	//

	BWF_Fitting_Results MorstinLinearDyBestWeightedResults;
	MorstinBWF.SetValues(std::vector<double> {7.87804976156578e-07, -2.13716006297374e-05, 0.000114200443445769, 5424.35985122075});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.0893680456034658});
	MorstinLinearDyBestWeightedResults.alphaFunc = MorstinBWF;
	MorstinLinearDyBestWeightedResults.betaFunc = FixedBWF;

	BWF_Fitting_Results MorstinLinearDyBestResults;
	MorstinBWF.SetValues(std::vector<double> {1.33638439092475e-06, -2.98606802538980e-05, 0.000137006580329845, 3010.90150042332});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.111377434741229});
	MorstinLinearDyBestResults.alphaFunc = MorstinBWF;
	MorstinLinearDyBestResults.betaFunc = FixedBWF;

	BWF_Fitting_Results CubicLinearDyBestRegularResults;
	// CubicBWF.SetValues(std::vector<double> {0.00838701932847768,0.185289234216632,-0.0250980413297348,0.000372587152280048});
	CubicBWF.SetValues(std::vector<double> {0.000337133584729527, -0.0222918654682546, 0.157099764860757, 0.0431408242731182});
	//FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.112398569166481});
	CubicLinearDyBestRegularResults.alphaFunc = CubicBWF;
	CubicLinearDyBestRegularResults.betaFunc = FixedBWF;

	BWF_Fitting_Results CubicLinearDyBestWeightedResults;
	// CubicBWF.SetValues(std::vector<double> {0.00838701932847768,0.185289234216632,-0.0250980413297348,0.000372587152280048});
	CubicBWF.SetValues(std::vector<double> {0.000372587152280048, -0.0250980413297348, 0.185289234216632, 0.00838701932847768});
	//FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.115994363406138});
	CubicLinearDyBestWeightedResults.alphaFunc = CubicBWF;
	CubicLinearDyBestWeightedResults.betaFunc = FixedBWF;

	BWF_Fitting_Results FourthLinearDyBestRegularResults;
	// CubicBWF.SetValues(std::vector<double> {0.00838701932847768,0.185289234216632,-0.0250980413297348,0.000372587152280048});
	FourthBWF.SetValues(std::vector<double> {1.55166949492021e-05, -0.00133016655165665, 0.0108472048475535, 1.32416754255453e-07, 0.171605600583206});
	//FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.114949030584930});
	FourthLinearDyBestRegularResults.alphaFunc = FourthBWF;
	FourthLinearDyBestRegularResults.betaFunc = FixedBWF;

	BWF_Fitting_Results FourthLETBest;
	FourthBWF.SetValues(std::vector<double> {1.14779e-05, -1.47913e-05, 0.000879819, -0.00200952, 0.212864});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.114996});
	FourthLETBest.alphaFunc = FourthBWF;
	FourthLETBest.betaFunc = FixedBWF;


	BWF_Fitting_Results CubicLETBest;
	CubicBWF.SetValues(std::vector<double> {0.000405226, -0.00394187, 0.0168077, 0.193427});
	FixedBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return params[0];}, 1);
	FixedBWF.SetValues(std::vector<double> {0.114966});
	CubicLETBest.alphaFunc = CubicBWF;
	CubicLETBest.betaFunc = FixedBWF;

	QEBWF.SetValues(std::vector<double>  {  0.209024, -0.105069, 0.000636537 });
	FixedBWF.SetValues(std::vector<double> {0.114731});
	BWF_Fitting_Results BestLETH460;
	BestLETH460.alphaFunc = QEBWF;
	BestLETH460.betaFunc = FixedBWF;

	AlphaBeta_Fitter fitter{};
	fitter.SetCellStudyParameters(H460FittingParams);
	double* AlphaBeta = fitter.Fit(nullptr,false);

	SurvivalDataMultigraph(c, legend, H460FittingParams);
	// TAttLine lineStyle{};
	lineStyle.SetLineColor(kRed+2);
	lineStyle.SetLineWidth(5);
	lineStyle.SetLineStyle(1);
	GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "Fourth - d(y)", H460FittingParams, FourthLinearDyBestRegularResults);
	lineStyle.SetLineColor(kBlue+2);
	GeneralizedBWFMultigraphPlotterLETAlphaBeta(c, legend, lineStyle, "Fourth - LET", H460FittingParams, FourthLETBest);
	lineStyle.SetLineColor(kGreen+3);
	MultigraphSurvivalFitPlotter(c, legend, lineStyle, "Direct Fit", H460FittingParams, AlphaBeta);
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/DirectFitH460.jpg";
	// GeneralizedBWFMultigraphPlotterLETAlphaBeta(c, legend, lineStyle, "QE - LET", H460FittingParams, BestLETH460);
	// GeneralizedBWFMultigraphPlotterAlphaBeta(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/survival_H460_Aug23_FourthLETvsDYvsDirect.jpg";
	c->SaveAs((TString)outputName); 

}

void H460_Consecutive_Fitting()
{
	SetupH460SurvivalParameters();

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
	// Fixed Beta Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << std::endl << "Linear, Fixed" << std::endl;
	// LinearBWF.SetValues(std::vector<double> {0.0399, 0.0457});
	LinearBWF.SetValues(std::vector<double> {0.0457, 0.0399});
	FixedBWF.SetValues(std::vector<double> {0.1355}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	double penaltyWeight = 0.0; 
	BWFFitter.SetPositiveConstrained(true, penaltyWeight);
	BWFFitter.Initialize();
	auto linearFixedResults = BWFFitter.Fit();
	linearFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFixedResults.betaFunc);

	std::cout << std::endl << "Quadratic, Fixed" << std::endl;
	std::vector<double> quadraticGuesses;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.push_back(0);
	quadraticGuesses.push_back(linearFixedResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(linearFixedResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticFixedResults = BWFFitter.Fit();
	quadraticFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFixedResults.betaFunc);

	std::cout << std::endl << "Cubic, Fixed" << std::endl;
	std::vector<double> cubicGuesses;
	cubicGuesses.push_back(0);
	cubicGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicFixedResults = BWFFitter.Fit();
	cubicFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFixedResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial, Fixed" << std::endl;
	std::vector<double> fourthGuesses;
	fourthGuesses.push_back(0);
	fourthGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthFixedResults = BWFFitter.Fit();
	fourthFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFixedResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial, Fixed" << std::endl;
	std::vector<double> fifthGuesses;
	fifthGuesses.push_back(0);
	fifthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthFixedResults = BWFFitter.Fit();
	fifthFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFixedResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q, Fixed" << std::endl;
	std::vector<double> QGuess;
	QGuess.push_back(linearFixedResults.alphaFunc.GetFittingParams()[2]);
	QGuess.push_back(0);
	QBWF.SetValues(std::vector<double> {0, 0.0001});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QFixedResults = BWFFitter.Fit();
	QFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFixedResults.betaFunc);

	std::cout << std::endl << "LE, Fixed" << std::endl;
	std::vector<double> LEGuesses;
	//Set alpha guess
	LEGuesses.push_back(0);
	LEGuesses.push_back(0.1);
	LEGuesses.push_back(linearFixedResults.alphaFunc.GetFittingParams()[1]);
	LEBWF.SetValues(std::vector<double> {LEGuesses});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	//Set beta guess
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LEFixedResults = BWFFitter.Fit();
	LEFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFixedResults.betaFunc);

	std::cout << std::endl << "QE, Fixed" << std::endl;
	std::vector<double> QEGuesses;
	QEGuesses.push_back(0);
	QEGuesses.push_back(0.1);
	QEGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[2]);
	QEBWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QEFixedResults = BWFFitter.Fit();
	QEFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFixedResults.betaFunc);

	std::cout << std::endl << "LQE, Fixed" << std::endl;
	std::vector<double> LQEGuesses;
	LQEGuesses.push_back(0);
	LQEGuesses.push_back(0.1);
	LQEGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[2]);
	LQEGuesses.push_back(linearFixedResults.alphaFunc.GetFittingParams()[1]);
	LQEBWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize(); 
	auto LQEFixedResults = BWFFitter.Fit();
	LQEFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFixedResults.betaFunc);

	std::cout << std::endl << "LE2, Fixed" << std::endl;
	LE2BWF.SetValues(LEGuesses); //same starting point as LE
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LE2FixedResults = BWFFitter.Fit();
	LE2FixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FixedResults.betaFunc);

	std::cout << std::endl << "QE2, Fixed" << std::endl;
	QE2BWF.SetValues(QEGuesses);
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto QE2FixedResults = BWFFitter.Fit();
	QE2FixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FixedResults.betaFunc);

	std::cout << std::endl << "LQE2, Fixed" << std::endl;
	LQE2BWF.SetValues(LQEGuesses);
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto LQE2FixedResults = BWFFitter.Fit();
	LQE2FixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FixedResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin, Fixed" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {1.73008e-07, -5.53241e-06, 2.53421e-05, 13541.4}); //{2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto MorstinFixedResults = BWFFitter.Fit();
	MorstinFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFixedResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset, Fixed" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 40, 80, 3});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// FixedBWF.SetValues(std::vector<double> {0});
	// BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetFixedResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetFixedResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFixedResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFixedResults.betaFunc);

	//
	//Linear Beta Fitting
	//

		//
		// The Polynomials
		//

	std::cout << std::endl << "Linear, Linear" << std::endl;
	LinearBWF.SetValues(std::vector<double> {linearFixedResults.alphaFunc.GetFittingParams()[1], linearFixedResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	LinearBWF.SetValues(std::vector<double> {0, linearFixedResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearLinearResults = BWFFitter.Fit();
	linearLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearLinearResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.clear();
	quadraticGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[2]);
	quadraticGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(quadraticFixedResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticLinearResults = BWFFitter.Fit();
	quadraticLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticLinearResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	cubicGuesses.clear();
	cubicGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[3]);
	cubicGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(cubicFixedResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicLinearResults = BWFFitter.Fit();
	cubicLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicLinearResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	fourthGuesses.clear();
	fourthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[4]);
	fourthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(fourthFixedResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthLinearResults = BWFFitter.Fit();
	fourthLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthLinearResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	fifthGuesses.clear();
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[5]);
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fifthFixedResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthLinearResults = BWFFitter.Fit();
	fifthLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthLinearResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	QBWF.SetValues(std::vector<double>{QFixedResults.alphaFunc.GetFittingParams()[1],QFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	LinearBWF.SetValues(std::vector<double> {0, QFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QLinearResults = BWFFitter.Fit();
	QLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QLinearResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {LEFixedResults.alphaFunc.GetFittingParams()[2],LEFixedResults.alphaFunc.GetFittingParams()[1],LEFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	LinearBWF.SetValues(std::vector<double> {0, LEFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LELinearResults = BWFFitter.Fit();
	LELinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LELinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LELinearResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {QEFixedResults.alphaFunc.GetFittingParams()[2],QEFixedResults.alphaFunc.GetFittingParams()[1],QEFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	LinearBWF.SetValues(std::vector<double> {0,QEFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QELinearResults = BWFFitter.Fit();
	QELinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QELinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QELinearResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	LQEBWF.SetValues(std::vector<double> {LQEFixedResults.alphaFunc.GetFittingParams()[3],LQEFixedResults.alphaFunc.GetFittingParams()[2],LQEFixedResults.alphaFunc.GetFittingParams()[1],LQEFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	LinearBWF.SetValues(std::vector<double> {0,LQEFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize(); 
	auto LQELinearResults = BWFFitter.Fit();
	LQELinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQELinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQELinearResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double> {LE2FixedResults.alphaFunc.GetFittingParams()[2],LE2FixedResults.alphaFunc.GetFittingParams()[1],LE2FixedResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	LinearBWF.SetValues(std::vector<double> {0, LE2FixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LE2LinearResults = BWFFitter.Fit();
	LE2LinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2LinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2LinearResults.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {QE2FixedResults.alphaFunc.GetFittingParams()[2],QE2FixedResults.alphaFunc.GetFittingParams()[1],QE2FixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	LinearBWF.SetValues(std::vector<double> {0, QE2FixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto QE2LinearResults = BWFFitter.Fit();
	QE2LinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2LinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2LinearResults.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {LQE2FixedResults.alphaFunc.GetFittingParams()[3],LQE2FixedResults.alphaFunc.GetFittingParams()[2],LQE2FixedResults.alphaFunc.GetFittingParams()[1],LQE2FixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	LinearBWF.SetValues(std::vector<double> {0, LQE2FixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto LQE2LinearResults = BWFFitter.Fit();
	LQE2LinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2LinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2LinearResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {MorstinFixedResults.alphaFunc.GetFittingParams()[3],MorstinFixedResults.alphaFunc.GetFittingParams()[2],MorstinFixedResults.alphaFunc.GetFittingParams()[1],MorstinFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	LinearBWF.SetValues(std::vector<double> {0, MorstinFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto MorstinLinearResults = BWFFitter.Fit();
	MorstinLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinLinearResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {GaussianVarAmplitudePlusOffsetFixedResults.alphaFunc.GetFittingParams()[3],GaussianVarAmplitudePlusOffsetFixedResults.alphaFunc.GetFittingParams()[2],GaussianVarAmplitudePlusOffsetFixedResults.alphaFunc.GetFittingParams()[1],GaussianVarAmplitudePlusOffsetFixedResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	LinearBWF.SetValues(std::vector<double> {0, GaussianVarAmplitudePlusOffsetFixedResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetLinearResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetLinearResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetLinearResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetLinearResults.betaFunc);

	//
	//Quadratic Beta Fitting
	//

		//
		// The Polynomials
		//

	std::cout << std::endl << "Linear, Quadratic" << std::endl;
	LinearBWF.SetValues(std::vector<double> {linearLinearResults.alphaFunc.GetFittingParams()[1], linearLinearResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	QuadraticBWF.SetValues(std::vector<double> {0, linearLinearResults.betaFunc.GetFittingParams()[1], linearLinearResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearQuadraticResults = BWFFitter.Fit();
	linearQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearQuadraticResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.clear();
	quadraticGuesses.push_back(quadraticLinearResults.alphaFunc.GetFittingParams()[2]);
	quadraticGuesses.push_back(quadraticLinearResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(quadraticLinearResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	// QuadraticBWF.SetValues(std::vector<double> {0, quadraticLinearResults.betaFunc.GetFittingParams()[1], quadraticLinearResults.betaFunc.GetFittingParams()[0]}); 
	// BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticQuadraticResults = BWFFitter.Fit();
	quadraticQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticQuadraticResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	cubicGuesses.clear();
	cubicGuesses.push_back(cubicLinearResults.alphaFunc.GetFittingParams()[3]);
	cubicGuesses.push_back(cubicLinearResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(cubicLinearResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(cubicLinearResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	// QuadraticBWF.SetValues(std::vector<double> {0, cubicLinearResults.betaFunc.GetFittingParams()[1], cubicLinearResults.betaFunc.GetFittingParams()[0]}); 
	// BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto cubicQuadraticResults = BWFFitter.Fit();
	cubicQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicQuadraticResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	fourthGuesses.clear();
	fourthGuesses.push_back(fourthLinearResults.alphaFunc.GetFittingParams()[4]);
	fourthGuesses.push_back(fourthLinearResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(fourthLinearResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(fourthLinearResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(fourthLinearResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	// QuadraticBWF.SetValues(std::vector<double> {0, fourthLinearResults.betaFunc.GetFittingParams()[1], fourthLinearResults.betaFunc.GetFittingParams()[0]}); 
	// BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto fourthQuadraticResults = BWFFitter.Fit();
	fourthQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthQuadraticResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	fifthGuesses.clear();
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[5]);
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fifthLinearResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	// QuadraticBWF.SetValues(std::vector<double> {0, fifthLinearResults.betaFunc.GetFittingParams()[1], fifthLinearResults.betaFunc.GetFittingParams()[0]}); 
	// BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto fifthQuadraticResults = BWFFitter.Fit();
	fifthQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthQuadraticResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	QBWF.SetValues(std::vector<double>{QLinearResults.alphaFunc.GetFittingParams()[1],QLinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	QuadraticBWF.SetValues(std::vector<double> {0, QLinearResults.betaFunc.GetFittingParams()[1], QLinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QQuadraticResults = BWFFitter.Fit();
	QQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QQuadraticResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {LELinearResults.alphaFunc.GetFittingParams()[2],LELinearResults.alphaFunc.GetFittingParams()[1],LELinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0, LELinearResults.betaFunc.GetFittingParams()[1], LELinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LEQuadraticResults = BWFFitter.Fit();
	LEQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEQuadraticResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {QELinearResults.alphaFunc.GetFittingParams()[2],QELinearResults.alphaFunc.GetFittingParams()[1],QELinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,  QELinearResults.betaFunc.GetFittingParams()[1], QELinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QEQuadraticResults = BWFFitter.Fit();
	QEQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEQuadraticResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	LQEBWF.SetValues(std::vector<double> {LQELinearResults.alphaFunc.GetFittingParams()[3],LQELinearResults.alphaFunc.GetFittingParams()[2],LQELinearResults.alphaFunc.GetFittingParams()[1],LQELinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	QuadraticBWF.SetValues(std::vector<double> {0,LQELinearResults.betaFunc.GetFittingParams()[1], LQELinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize(); 
	auto LQEQuadraticResults = BWFFitter.Fit();
	LQEQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEQuadraticResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double> {LE2LinearResults.alphaFunc.GetFittingParams()[2],LE2LinearResults.alphaFunc.GetFittingParams()[1],LE2LinearResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0, LE2LinearResults.betaFunc.GetFittingParams()[1], LE2LinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LE2QuadraticResults = BWFFitter.Fit();
	LE2QuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2QuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2QuadraticResults.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {QE2LinearResults.alphaFunc.GetFittingParams()[2],QE2LinearResults.alphaFunc.GetFittingParams()[1],QE2LinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0, QE2LinearResults.betaFunc.GetFittingParams()[1], QE2LinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto QE2QuadraticResults = BWFFitter.Fit();
	QE2QuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2QuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2QuadraticResults.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {LQE2LinearResults.alphaFunc.GetFittingParams()[3],LQE2LinearResults.alphaFunc.GetFittingParams()[2],LQE2LinearResults.alphaFunc.GetFittingParams()[1],LQE2LinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	QuadraticBWF.SetValues(std::vector<double> {0, LQE2LinearResults.betaFunc.GetFittingParams()[1], LQE2LinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto LQE2QuadraticResults = BWFFitter.Fit();
	LQE2QuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2QuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2QuadraticResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {MorstinLinearResults.alphaFunc.GetFittingParams()[3],MorstinLinearResults.alphaFunc.GetFittingParams()[2],MorstinLinearResults.alphaFunc.GetFittingParams()[1],MorstinLinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	QuadraticBWF.SetValues(std::vector<double> {0, MorstinLinearResults.betaFunc.GetFittingParams()[1], MorstinLinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto MorstinQuadraticResults = BWFFitter.Fit();
	MorstinQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinQuadraticResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {GaussianVarAmplitudePlusOffsetLinearResults.alphaFunc.GetFittingParams()[3],GaussianVarAmplitudePlusOffsetLinearResults.alphaFunc.GetFittingParams()[2],GaussianVarAmplitudePlusOffsetLinearResults.alphaFunc.GetFittingParams()[1],GaussianVarAmplitudePlusOffsetLinearResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	QuadraticBWF.SetValues(std::vector<double> {0, GaussianVarAmplitudePlusOffsetLinearResults.betaFunc.GetFittingParams()[1], GaussianVarAmplitudePlusOffsetLinearResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetQuadraticResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetQuadraticResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetQuadraticResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetQuadraticResults.betaFunc);

	//
	//Cubic Beta Fitting
	//

		//
		// The Polynomials
		//

	std::cout << std::endl << "Linear, Cubic" << std::endl;
	LinearBWF.SetValues(std::vector<double> {linearQuadraticResults.alphaFunc.GetFittingParams()[1], linearQuadraticResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	CubicBWF.SetValues(std::vector<double> {0, linearQuadraticResults.betaFunc.GetFittingParams()[2], linearQuadraticResults.betaFunc.GetFittingParams()[1], linearQuadraticResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearCubicResults = BWFFitter.Fit();
	linearCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearCubicResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.clear();
	quadraticGuesses.push_back(quadraticQuadraticResults.alphaFunc.GetFittingParams()[2]);
	quadraticGuesses.push_back(quadraticQuadraticResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(quadraticQuadraticResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	CubicBWF.SetValues(std::vector<double> {0, quadraticQuadraticResults.betaFunc.GetFittingParams()[2], quadraticQuadraticResults.betaFunc.GetFittingParams()[1], quadraticQuadraticResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto quadraticCubicResults = BWFFitter.Fit();
	quadraticCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticCubicResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	cubicGuesses.clear();
	cubicGuesses.push_back(cubicQuadraticResults.alphaFunc.GetFittingParams()[3]);
	cubicGuesses.push_back(cubicQuadraticResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(cubicQuadraticResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(cubicQuadraticResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicCubicResults = BWFFitter.Fit();
	cubicCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicCubicResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	fourthGuesses.clear();
	fourthGuesses.push_back(fourthQuadraticResults.alphaFunc.GetFittingParams()[4]);
	fourthGuesses.push_back(fourthQuadraticResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(fourthQuadraticResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(fourthQuadraticResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(fourthQuadraticResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthCubicResults = BWFFitter.Fit();
	fourthCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthCubicResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	fifthGuesses.clear();
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[5]);
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fifthQuadraticResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthCubicResults = BWFFitter.Fit();
	fifthCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthCubicResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	QBWF.SetValues(std::vector<double>{QQuadraticResults.alphaFunc.GetFittingParams()[1],QQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	CubicBWF.SetValues(std::vector<double> {0, QQuadraticResults.betaFunc.GetFittingParams()[2], QQuadraticResults.betaFunc.GetFittingParams()[1], QQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QCubicResults = BWFFitter.Fit();
	QCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QCubicResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {LEQuadraticResults.alphaFunc.GetFittingParams()[2],LEQuadraticResults.alphaFunc.GetFittingParams()[1],LEQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	CubicBWF.SetValues(std::vector<double> {0, LEQuadraticResults.betaFunc.GetFittingParams()[2], LEQuadraticResults.betaFunc.GetFittingParams()[1], LEQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LECubicResults = BWFFitter.Fit();
	LECubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LECubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LECubicResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {QEQuadraticResults.alphaFunc.GetFittingParams()[2],QEQuadraticResults.alphaFunc.GetFittingParams()[1],QEQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	CubicBWF.SetValues(std::vector<double> {0,  QEQuadraticResults.betaFunc.GetFittingParams()[2], QEQuadraticResults.betaFunc.GetFittingParams()[1], QEQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QECubicResults = BWFFitter.Fit();
	QECubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QECubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QECubicResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	LQEBWF.SetValues(std::vector<double> {LQEQuadraticResults.alphaFunc.GetFittingParams()[3],LQEQuadraticResults.alphaFunc.GetFittingParams()[2],LQEQuadraticResults.alphaFunc.GetFittingParams()[1],LQEQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	CubicBWF.SetValues(std::vector<double> {0,LQEQuadraticResults.betaFunc.GetFittingParams()[2], LQEQuadraticResults.betaFunc.GetFittingParams()[1], LQEQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize(); 
	auto LQECubicResults = BWFFitter.Fit();
	LQECubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQECubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQECubicResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double> {LE2QuadraticResults.alphaFunc.GetFittingParams()[2],LE2QuadraticResults.alphaFunc.GetFittingParams()[1],LE2QuadraticResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	CubicBWF.SetValues(std::vector<double> {0, LE2QuadraticResults.betaFunc.GetFittingParams()[2], LE2QuadraticResults.betaFunc.GetFittingParams()[1], LE2QuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LE2CubicResults = BWFFitter.Fit();
	LE2CubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2CubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2CubicResults.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {QE2QuadraticResults.alphaFunc.GetFittingParams()[2],QE2QuadraticResults.alphaFunc.GetFittingParams()[1],QE2QuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	CubicBWF.SetValues(std::vector<double> {0, QE2QuadraticResults.betaFunc.GetFittingParams()[2], QE2QuadraticResults.betaFunc.GetFittingParams()[1], QE2QuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto QE2CubicResults = BWFFitter.Fit();
	QE2CubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2CubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2CubicResults.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {LQE2QuadraticResults.alphaFunc.GetFittingParams()[3],LQE2QuadraticResults.alphaFunc.GetFittingParams()[2],LQE2QuadraticResults.alphaFunc.GetFittingParams()[1],LQE2QuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	CubicBWF.SetValues(std::vector<double> {0, LQE2QuadraticResults.betaFunc.GetFittingParams()[2], LQE2QuadraticResults.betaFunc.GetFittingParams()[1], LQE2QuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto LQE2CubicResults = BWFFitter.Fit();
	LQE2CubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2CubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2CubicResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {MorstinQuadraticResults.alphaFunc.GetFittingParams()[3],MorstinQuadraticResults.alphaFunc.GetFittingParams()[2],MorstinQuadraticResults.alphaFunc.GetFittingParams()[1],MorstinQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	CubicBWF.SetValues(std::vector<double> {0, MorstinQuadraticResults.betaFunc.GetFittingParams()[2], MorstinQuadraticResults.betaFunc.GetFittingParams()[1], MorstinQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto MorstinCubicResults = BWFFitter.Fit();
	MorstinCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinCubicResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {GaussianVarAmplitudePlusOffsetQuadraticResults.alphaFunc.GetFittingParams()[3],GaussianVarAmplitudePlusOffsetQuadraticResults.alphaFunc.GetFittingParams()[2],GaussianVarAmplitudePlusOffsetQuadraticResults.alphaFunc.GetFittingParams()[1],GaussianVarAmplitudePlusOffsetQuadraticResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	CubicBWF.SetValues(std::vector<double> {0, GaussianVarAmplitudePlusOffsetQuadraticResults.betaFunc.GetFittingParams()[2], GaussianVarAmplitudePlusOffsetQuadraticResults.betaFunc.GetFittingParams()[1], GaussianVarAmplitudePlusOffsetQuadraticResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetCubicResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetCubicResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetCubicResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetCubicResults.betaFunc);

	//
	//Fourth Beta Fitting
	//

		//
		// The Polynomials
		//

	std::cout << std::endl << "Linear, Fourth" << std::endl;
	LinearBWF.SetValues(std::vector<double> {linearCubicResults.alphaFunc.GetFittingParams()[1], linearCubicResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FourthBWF.SetValues(std::vector<double> {0, linearCubicResults.betaFunc.GetFittingParams()[3], linearCubicResults.betaFunc.GetFittingParams()[2], linearCubicResults.betaFunc.GetFittingParams()[1], linearCubicResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearFourthResults = BWFFitter.Fit();
	linearFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFourthResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.clear();
	quadraticGuesses.push_back(quadraticCubicResults.alphaFunc.GetFittingParams()[2]);
	quadraticGuesses.push_back(quadraticCubicResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(quadraticCubicResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	FourthBWF.SetValues(std::vector<double> {0, quadraticCubicResults.betaFunc.GetFittingParams()[3], quadraticCubicResults.betaFunc.GetFittingParams()[2], quadraticCubicResults.betaFunc.GetFittingParams()[1], quadraticCubicResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	auto quadraticFourthResults = BWFFitter.Fit();
	quadraticFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFourthResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	cubicGuesses.clear();
	cubicGuesses.push_back(cubicCubicResults.alphaFunc.GetFittingParams()[3]);
	cubicGuesses.push_back(cubicCubicResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(cubicCubicResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(cubicCubicResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicFourthResults = BWFFitter.Fit();
	cubicFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFourthResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	fourthGuesses.clear();
	fourthGuesses.push_back(fourthCubicResults.alphaFunc.GetFittingParams()[4]);
	fourthGuesses.push_back(fourthCubicResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(fourthCubicResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(fourthCubicResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(fourthCubicResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthFourthResults = BWFFitter.Fit();
	fourthFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFourthResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	fifthGuesses.clear();
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[5]);
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fifthCubicResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthFourthResults = BWFFitter.Fit();
	fifthFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFourthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	QBWF.SetValues(std::vector<double>{QCubicResults.alphaFunc.GetFittingParams()[1],QCubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FourthBWF.SetValues(std::vector<double> {0, QCubicResults.betaFunc.GetFittingParams()[3], QCubicResults.betaFunc.GetFittingParams()[2], QCubicResults.betaFunc.GetFittingParams()[1], QCubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QFourthResults = BWFFitter.Fit();
	QFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFourthResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {LECubicResults.alphaFunc.GetFittingParams()[2],LECubicResults.alphaFunc.GetFittingParams()[1],LECubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	FourthBWF.SetValues(std::vector<double> {0, LECubicResults.betaFunc.GetFittingParams()[3], LECubicResults.betaFunc.GetFittingParams()[2], LECubicResults.betaFunc.GetFittingParams()[1], LECubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LEFourthResults = BWFFitter.Fit();
	LEFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFourthResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {QECubicResults.alphaFunc.GetFittingParams()[2],QECubicResults.alphaFunc.GetFittingParams()[1],QECubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	FourthBWF.SetValues(std::vector<double> {0,  QECubicResults.betaFunc.GetFittingParams()[3], QECubicResults.betaFunc.GetFittingParams()[2], QECubicResults.betaFunc.GetFittingParams()[1], QECubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QEFourthResults = BWFFitter.Fit();
	QEFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFourthResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	LQEBWF.SetValues(std::vector<double> {LQECubicResults.alphaFunc.GetFittingParams()[3],LQECubicResults.alphaFunc.GetFittingParams()[2],LQECubicResults.alphaFunc.GetFittingParams()[1],LQECubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	FourthBWF.SetValues(std::vector<double> {0, LQECubicResults.betaFunc.GetFittingParams()[3], LQECubicResults.betaFunc.GetFittingParams()[2], LQECubicResults.betaFunc.GetFittingParams()[1], LQECubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize(); 
	auto LQEFourthResults = BWFFitter.Fit();
	LQEFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFourthResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double> {LE2CubicResults.alphaFunc.GetFittingParams()[2],LE2CubicResults.alphaFunc.GetFittingParams()[1],LE2CubicResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FourthBWF.SetValues(std::vector<double> {0, LE2CubicResults.betaFunc.GetFittingParams()[3], LE2CubicResults.betaFunc.GetFittingParams()[2], LE2CubicResults.betaFunc.GetFittingParams()[1], LE2CubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LE2FourthResults = BWFFitter.Fit();
	LE2FourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FourthResults.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {QE2CubicResults.alphaFunc.GetFittingParams()[2],QE2CubicResults.alphaFunc.GetFittingParams()[1],QE2CubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FourthBWF.SetValues(std::vector<double> {0, QE2CubicResults.betaFunc.GetFittingParams()[3], QE2CubicResults.betaFunc.GetFittingParams()[2], QE2CubicResults.betaFunc.GetFittingParams()[1], QE2CubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto QE2FourthResults = BWFFitter.Fit();
	QE2FourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FourthResults.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {LQE2CubicResults.alphaFunc.GetFittingParams()[3],LQE2CubicResults.alphaFunc.GetFittingParams()[2],LQE2CubicResults.alphaFunc.GetFittingParams()[1],LQE2CubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FourthBWF.SetValues(std::vector<double> {0, LQE2CubicResults.betaFunc.GetFittingParams()[3], LQE2CubicResults.betaFunc.GetFittingParams()[2], LQE2CubicResults.betaFunc.GetFittingParams()[1], LQE2CubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto LQE2FourthResults = BWFFitter.Fit();
	LQE2FourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FourthResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {MorstinCubicResults.alphaFunc.GetFittingParams()[3],MorstinCubicResults.alphaFunc.GetFittingParams()[2],MorstinCubicResults.alphaFunc.GetFittingParams()[1],MorstinCubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FourthBWF.SetValues(std::vector<double> {0, MorstinCubicResults.betaFunc.GetFittingParams()[3], MorstinCubicResults.betaFunc.GetFittingParams()[2], MorstinCubicResults.betaFunc.GetFittingParams()[1], MorstinCubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto MorstinFourthResults = BWFFitter.Fit();
	MorstinFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFourthResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {GaussianVarAmplitudePlusOffsetCubicResults.alphaFunc.GetFittingParams()[3],GaussianVarAmplitudePlusOffsetCubicResults.alphaFunc.GetFittingParams()[2],GaussianVarAmplitudePlusOffsetCubicResults.alphaFunc.GetFittingParams()[1],GaussianVarAmplitudePlusOffsetCubicResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FourthBWF.SetValues(std::vector<double> {0, GaussianVarAmplitudePlusOffsetCubicResults.betaFunc.GetFittingParams()[3], GaussianVarAmplitudePlusOffsetCubicResults.betaFunc.GetFittingParams()[2], GaussianVarAmplitudePlusOffsetCubicResults.betaFunc.GetFittingParams()[1], GaussianVarAmplitudePlusOffsetCubicResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetFourthResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetFourthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFourthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFourthResults.betaFunc);

	//
	//Fifth Beta Fitting
	//

		//
		// The Polynomials
		//

	std::cout << std::endl << "Linear, Fifth" << std::endl;
	LinearBWF.SetValues(std::vector<double> {linearFourthResults.alphaFunc.GetFittingParams()[1], linearFourthResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FifthBWF.SetValues(std::vector<double> {0, linearFourthResults.betaFunc.GetFittingParams()[4], linearFourthResults.betaFunc.GetFittingParams()[3], linearFourthResults.betaFunc.GetFittingParams()[2], linearFourthResults.betaFunc.GetFittingParams()[1], linearFourthResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearFifthResults = BWFFitter.Fit();
	linearFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(linearFifthResults.betaFunc);

	std::cout << std::endl << "Quadratic" << std::endl;
	//Take the results from the last fit and use them to inform this fit
	quadraticGuesses.clear();
	quadraticGuesses.push_back(quadraticFourthResults.alphaFunc.GetFittingParams()[2]);
	quadraticGuesses.push_back(quadraticFourthResults.alphaFunc.GetFittingParams()[1]);
	quadraticGuesses.push_back(quadraticFourthResults.alphaFunc.GetFittingParams()[0]);
	QuadraticBWF.SetValues(quadraticGuesses);
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	FifthBWF.SetValues(std::vector<double> {0, quadraticFourthResults.betaFunc.GetFittingParams()[4], quadraticFourthResults.betaFunc.GetFittingParams()[3], quadraticFourthResults.betaFunc.GetFittingParams()[2], quadraticFourthResults.betaFunc.GetFittingParams()[1], quadraticFourthResults.betaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto quadraticFifthResults = BWFFitter.Fit();
	quadraticFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(quadraticFifthResults.betaFunc);

	std::cout << std::endl << "Cubic" << std::endl;
	cubicGuesses.clear();
	cubicGuesses.push_back(cubicFourthResults.alphaFunc.GetFittingParams()[3]);
	cubicGuesses.push_back(cubicFourthResults.alphaFunc.GetFittingParams()[2]);
	cubicGuesses.push_back(cubicFourthResults.alphaFunc.GetFittingParams()[1]);
	cubicGuesses.push_back(cubicFourthResults.alphaFunc.GetFittingParams()[0]);
	CubicBWF.SetValues(cubicGuesses);
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	auto cubicFifthResults = BWFFitter.Fit();
	cubicFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(cubicFifthResults.betaFunc);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	fourthGuesses.clear();
	fourthGuesses.push_back(fourthFourthResults.alphaFunc.GetFittingParams()[4]);
	fourthGuesses.push_back(fourthFourthResults.alphaFunc.GetFittingParams()[3]);
	fourthGuesses.push_back(fourthFourthResults.alphaFunc.GetFittingParams()[2]);
	fourthGuesses.push_back(fourthFourthResults.alphaFunc.GetFittingParams()[1]);
	fourthGuesses.push_back(fourthFourthResults.alphaFunc.GetFittingParams()[0]);
	FourthBWF.SetValues(fourthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	auto fourthFifthResults = BWFFitter.Fit();
	fourthFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fourthFifthResults.betaFunc);

	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	fifthGuesses.clear();
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[5]);
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[4]);
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[3]);
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[2]);
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[1]);
	fifthGuesses.push_back(fifthFourthResults.alphaFunc.GetFittingParams()[0]);
	FifthBWF.SetValues(fifthGuesses);
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthFifthResults = BWFFitter.Fit();
	fifthFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(fifthFifthResults.betaFunc);

		//
		// Mairani et al.
		//

	std::cout << std::endl << "Q" << std::endl;
	QBWF.SetValues(std::vector<double>{QFourthResults.alphaFunc.GetFittingParams()[1],QFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FifthBWF.SetValues(std::vector<double> {0, QFourthResults.betaFunc.GetFittingParams()[4], QFourthResults.betaFunc.GetFittingParams()[3], QFourthResults.betaFunc.GetFittingParams()[2], QFourthResults.betaFunc.GetFittingParams()[1], QFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QFifthResults = BWFFitter.Fit();
	QFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QFifthResults.betaFunc);

	std::cout << std::endl << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {LEFourthResults.alphaFunc.GetFittingParams()[2],LEFourthResults.alphaFunc.GetFittingParams()[1],LEFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	FifthBWF.SetValues(std::vector<double> {0, LEFourthResults.betaFunc.GetFittingParams()[4], LEFourthResults.betaFunc.GetFittingParams()[3], LEFourthResults.betaFunc.GetFittingParams()[2], LEFourthResults.betaFunc.GetFittingParams()[1], LEFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LEFifthResults = BWFFitter.Fit();
	LEFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LEFifthResults.betaFunc);

	std::cout << std::endl << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {QEFourthResults.alphaFunc.GetFittingParams()[2],QEFourthResults.alphaFunc.GetFittingParams()[1],QEFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	FifthBWF.SetValues(std::vector<double> {0,  QEFourthResults.betaFunc.GetFittingParams()[4], QEFourthResults.betaFunc.GetFittingParams()[3], QEFourthResults.betaFunc.GetFittingParams()[2], QEFourthResults.betaFunc.GetFittingParams()[1], QEFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QEFifthResults = BWFFitter.Fit();
	QEFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QEFifthResults.betaFunc);

	std::cout << std::endl << "LQE" << std::endl;
	LQEBWF.SetValues(std::vector<double> {LQEFourthResults.alphaFunc.GetFittingParams()[3],LQEFourthResults.alphaFunc.GetFittingParams()[2],LQEFourthResults.alphaFunc.GetFittingParams()[1],LQEFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQEBWF);
	FifthBWF.SetValues(std::vector<double> {0, LQEFourthResults.betaFunc.GetFittingParams()[4], LQEFourthResults.betaFunc.GetFittingParams()[3], LQEFourthResults.betaFunc.GetFittingParams()[2], LQEFourthResults.betaFunc.GetFittingParams()[1], LQEFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize(); 
	auto LQEFifthResults = BWFFitter.Fit();
	LQEFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQEFifthResults.betaFunc);

	std::cout << std::endl << "LE2" << std::endl;
	LE2BWF.SetValues(std::vector<double> {LE2FourthResults.alphaFunc.GetFittingParams()[2],LE2FourthResults.alphaFunc.GetFittingParams()[1],LE2FourthResults.alphaFunc.GetFittingParams()[0]}); 
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	FifthBWF.SetValues(std::vector<double> {0, LE2FourthResults.betaFunc.GetFittingParams()[4], LE2FourthResults.betaFunc.GetFittingParams()[3], LE2FourthResults.betaFunc.GetFittingParams()[2], LE2FourthResults.betaFunc.GetFittingParams()[1], LE2FourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LE2FifthResults = BWFFitter.Fit();
	LE2FifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LE2FifthResults.betaFunc);

	std::cout << std::endl << "QE2" << std::endl;
	QE2BWF.SetValues(std::vector<double> {QE2FourthResults.alphaFunc.GetFittingParams()[2],QE2FourthResults.alphaFunc.GetFittingParams()[1],QE2FourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	FifthBWF.SetValues(std::vector<double> {0, QE2FourthResults.betaFunc.GetFittingParams()[4], QE2FourthResults.betaFunc.GetFittingParams()[3], QE2FourthResults.betaFunc.GetFittingParams()[2], QE2FourthResults.betaFunc.GetFittingParams()[1], QE2FourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto QE2FifthResults = BWFFitter.Fit();
	QE2FifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(QE2FifthResults.betaFunc);

	std::cout << std::endl << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {LQE2FourthResults.alphaFunc.GetFittingParams()[3],LQE2FourthResults.alphaFunc.GetFittingParams()[2],LQE2FourthResults.alphaFunc.GetFittingParams()[1],LQE2FourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	FifthBWF.SetValues(std::vector<double> {0, LQE2FourthResults.betaFunc.GetFittingParams()[4], LQE2FourthResults.betaFunc.GetFittingParams()[3], LQE2FourthResults.betaFunc.GetFittingParams()[2], LQE2FourthResults.betaFunc.GetFittingParams()[1], LQE2FourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto LQE2FifthResults = BWFFitter.Fit();
	LQE2FifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(LQE2FifthResults.betaFunc);

		//
		//Fitting the other functional forms
		//

	std::cout << std::endl << "Morstin" << std::endl;
	MorstinBWF.SetValues(std::vector<double> {MorstinFourthResults.alphaFunc.GetFittingParams()[3],MorstinFourthResults.alphaFunc.GetFittingParams()[2],MorstinFourthResults.alphaFunc.GetFittingParams()[1],MorstinFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	FifthBWF.SetValues(std::vector<double> {0, MorstinFourthResults.betaFunc.GetFittingParams()[4], MorstinFourthResults.betaFunc.GetFittingParams()[3], MorstinFourthResults.betaFunc.GetFittingParams()[2], MorstinFourthResults.betaFunc.GetFittingParams()[1], MorstinFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto MorstinFifthResults = BWFFitter.Fit();
	MorstinFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(MorstinFifthResults.betaFunc);

	std::cout << std::endl << "Gaussian plus offset" << std::endl;
	GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {GaussianVarAmplitudePlusOffsetFourthResults.alphaFunc.GetFittingParams()[3],GaussianVarAmplitudePlusOffsetFourthResults.alphaFunc.GetFittingParams()[2],GaussianVarAmplitudePlusOffsetFourthResults.alphaFunc.GetFittingParams()[1],GaussianVarAmplitudePlusOffsetFourthResults.alphaFunc.GetFittingParams()[0]});
	BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	FifthBWF.SetValues(std::vector<double> {0, GaussianVarAmplitudePlusOffsetFourthResults.betaFunc.GetFittingParams()[4], GaussianVarAmplitudePlusOffsetFourthResults.betaFunc.GetFittingParams()[3], GaussianVarAmplitudePlusOffsetFourthResults.betaFunc.GetFittingParams()[2], GaussianVarAmplitudePlusOffsetFourthResults.betaFunc.GetFittingParams()[1], GaussianVarAmplitudePlusOffsetFourthResults.betaFunc.GetFittingParams()[0]});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto GaussianVarAmplitudePlusOffsetFifthResults = BWFFitter.Fit();
	GaussianVarAmplitudePlusOffsetFifthResults.PrintRMSEMinusPenaltyFunction(penaltyWeight, 83);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFifthResults.alphaFunc);
	Ceres_BWF_Fitter::CheckFunctionNegativity(GaussianVarAmplitudePlusOffsetFifthResults.betaFunc);

}

void H460_Ceres()
{
	Best_Model_Plotting();
	// Survival_Fitting();
}