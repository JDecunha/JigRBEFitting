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
#include "Ceres_BWF_Fitter.h"
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
	TFile* dySpectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraCellStudy_1e3nm.root");
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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


	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	BWFFitter.SetBetaWeightingFunction(FixedBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();

	std::cout << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0454, 0.0453});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();

	std::cout << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.00248214, -0.0583471, 0.275426});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();

	std::cout << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 0.000198751, -0.0122275, 0.0654629, 0.110813});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 1.86239e-05, -0.00178817, 0.0269139, -0.118944, 0.271133});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();

		//
		// Mairani et al.
		//

	std::cout << "Q" << std::endl;
	QBWF.SetValues(std::vector<double> {1});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QResults = BWFFitter.Fit();
	QResults.PrintBasic();

	std::cout << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {0.3, 0.1});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintBasic();

	std::cout << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {0.01, 1e-6});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintBasic();

	std::cout << "LQE" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	BWFFitter.Initialize(); //then initialize
	BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintBasic();

	std::cout << "LE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintBasic();

	std::cout << "QE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintBasic();

	std::cout << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {-0.000256818, -0.000171922, 0.0314693});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintBasic();

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

	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomial_CERES_BWFs_compto3rdpoly.jpg";
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
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_Ceres.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_Ceres.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_Ceres.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_Ceres.jpg";
	// c->SaveAs((TString)outputName); 


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

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Mairani_BWFs_constrained.jpg";
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
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_Ceres.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2.jpg";
	// c->SaveAs((TString)outputName); 

	// AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	// outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2.jpg";
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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


	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	BWFFitter.SetBetaWeightingFunction(LinearBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();

	std::cout << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0454, 0.0453});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();

	std::cout << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.00248214, -0.0583471, 0.275426});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();

	std::cout << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 0.000198751, -0.0122275, 0.0654629, 0.110813});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 1.86239e-05, -0.00178817, 0.0269139, -0.118944, 0.271133});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();

		//
		// Mairani et al.
		//

	std::cout << "Q" << std::endl;
	QBWF.SetValues(std::vector<double> {1});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QResults = BWFFitter.Fit();
	QResults.PrintBasic();

	std::cout << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {0.3, 0.1});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintBasic();

	std::cout << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {0.01, 1e-6});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintBasic();

	std::cout << "LQE" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	BWFFitter.Initialize(); //then initialize
	BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintBasic();

	std::cout << "LE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintBasic();

	std::cout << "QE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintBasic();

	std::cout << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {-0.000256818, -0.000171922, 0.0314693});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintBasic();

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

	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomial_CERES_BWFs_compto3rdpoly.jpg";
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
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_linear.jpg";
	c->SaveAs((TString)outputName); 


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

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Mairani_BWFs_constrained.jpg";
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

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_linear.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_linear.jpg";
	c->SaveAs((TString)outputName); 
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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


	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	BWFFitter.SetBetaWeightingFunction(QuadraticBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();

	std::cout << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0454, 0.0453});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();

	std::cout << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.00248214, -0.0583471, 0.275426});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();

	std::cout << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 0.000198751, -0.0122275, 0.0654629, 0.110813});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 1.86239e-05, -0.00178817, 0.0269139, -0.118944, 0.271133});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();

		//
		// Mairani et al.
		//

	std::cout << "Q" << std::endl;
	QBWF.SetValues(std::vector<double> {1});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QResults = BWFFitter.Fit();
	QResults.PrintBasic();

	std::cout << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {0.3, 0.1});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintBasic();

	std::cout << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {0.01, 1e-6});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintBasic();

	std::cout << "LQE" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	BWFFitter.Initialize(); //then initialize
	BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintBasic();

	std::cout << "LE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintBasic();

	std::cout << "QE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintBasic();

	std::cout << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {-0.000256818, -0.000171922, 0.0314693});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintBasic();

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

	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomial_CERES_BWFs_compto3rdpoly.jpg";
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
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_quadratic.jpg";
	c->SaveAs((TString)outputName); 


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

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Mairani_BWFs_constrained.jpg";
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

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_quadratic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_quadratic.jpg";
	c->SaveAs((TString)outputName); 
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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


	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	BWFFitter.SetBetaWeightingFunction(CubicBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();

	std::cout << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0454, 0.0453});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();

	std::cout << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.00248214, -0.0583471, 0.275426});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();

	std::cout << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 0.000198751, -0.0122275, 0.0654629, 0.110813});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 1.86239e-05, -0.00178817, 0.0269139, -0.118944, 0.271133});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();

		//
		// Mairani et al.
		//

	std::cout << "Q" << std::endl;
	QBWF.SetValues(std::vector<double> {1});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QResults = BWFFitter.Fit();
	QResults.PrintBasic();

	std::cout << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {0.3, 0.1});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintBasic();

	std::cout << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {0.01, 1e-6});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintBasic();

	std::cout << "LQE" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	BWFFitter.Initialize(); //then initialize
	BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintBasic();

	std::cout << "LE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintBasic();

	std::cout << "QE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintBasic();

	std::cout << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {-0.000256818, -0.000171922, 0.0314693});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintBasic();

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

	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomial_CERES_BWFs_compto3rdpoly.jpg";
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
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_cubic.jpg";
	c->SaveAs((TString)outputName); 


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

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Mairani_BWFs_constrained.jpg";
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

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_cubic.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_cubic.jpg";
	c->SaveAs((TString)outputName); 
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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


	//
	// The Actual Fitting
	//

		//
		// The Polynomials
		//

	Ceres_BWF_Fitter BWFFitter{};
	std::cout << "Linear" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LinearBWF);
	// FixedBWF.SetValues(std::vector<double> {0.13}); //Give beta value before passing into the fitter
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();

	std::cout << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0454, 0.0453});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();

	std::cout << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.00248214, -0.0583471, 0.275426});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();

	std::cout << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 0.000198751, -0.0122275, 0.0654629, 0.110813});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 1.86239e-05, -0.00178817, 0.0269139, -0.118944, 0.271133});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();

		//
		// Mairani et al.
		//

	std::cout << "Q" << std::endl;
	QBWF.SetValues(std::vector<double> {1});
	BWFFitter.SetAlphaWeightingFunction(QBWF);
	FourthBWF.SetValues(std::vector<double> {0.1,0.1,0.1,0.1,0.1});
	BWFFitter.SetBetaWeightingFunction(FourthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QResults = BWFFitter.Fit();
	QResults.PrintBasic();

	std::cout << "LE" << std::endl;
	LEBWF.SetValues(std::vector<double> {0.3, 0.1});
	BWFFitter.SetAlphaWeightingFunction(LEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LEResults = BWFFitter.Fit();
	LEResults.PrintBasic();

	std::cout << "QE" << std::endl;
	QEBWF.SetValues(std::vector<double> {0.01, 1e-6});
	BWFFitter.SetAlphaWeightingFunction(QEBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QEResults = BWFFitter.Fit();
	QEResults.PrintBasic();

	std::cout << "LQE" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	BWFFitter.Initialize(); //then initialize
	BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQEResults = BWFFitter.Fit();
	LQEResults.PrintBasic();

	std::cout << "LE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LE2Results = BWFFitter.Fit();
	LE2Results.PrintBasic();

	std::cout << "QE2" << std::endl;
	BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto QE2Results = BWFFitter.Fit();
	QE2Results.PrintBasic();

	std::cout << "LQE2" << std::endl;
	LQE2BWF.SetValues(std::vector<double> {-0.000256818, -0.000171922, 0.0314693});
	BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto LQE2Results = BWFFitter.Fit();
	LQE2Results.PrintBasic();

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

	// BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kBlue+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kRed+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	// lineStyle.SetLineColor(kOrange+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomial_CERES_BWFs_compto3rdpoly.jpg";
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
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_linear_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Quadratic", H460FittingParams, quadraticResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_quadratic_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Cubic", H460FittingParams, cubicResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_cubic_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fourth", H460FittingParams, fourthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fourth_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Fifth", H460FittingParams, fifthResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_fifth_fourth.jpg";
	c->SaveAs((TString)outputName); 


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

	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/Mairani_BWFs_constrained.jpg";
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

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "Q", H460FittingParams, QResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_Q_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE", H460FittingParams, LEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE", H460FittingParams, QEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE", H460FittingParams, LQEResults);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LE2", H460FittingParams, LE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LE2_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "QE2", H460FittingParams, QE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_QE2_fourth.jpg";
	c->SaveAs((TString)outputName); 

	AlphaBetaMultigraphResiduals(c, legend, lineStyle, "LQE2", H460FittingParams, LQE2Results);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/residuals_LQE2_fourth.jpg";
	c->SaveAs((TString)outputName); 
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
	QBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return (params[0]*linealEnergy*linealEnergy);}, 1);

	//Create LE BWF (linear-exponential)
	BiologicalWeightingFunction LEBWF;
	LEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);

	//Create QE BWF (linear-exponential)
	BiologicalWeightingFunction QEBWF;
	QEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy));}, 2);
	
	//Create LQE BWF (linear-exponential)
	BiologicalWeightingFunction LQEBWF;
	LQEBWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy));}, 3);
	
	//Create LE2 BWF
	BiologicalWeightingFunction LE2BWF;
	LE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create QE2 BWF
	BiologicalWeightingFunction QE2BWF;
	QE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy*linealEnergy)*std::exp(-params[1]*linealEnergy*linealEnergy));}, 2);
	
	//Create LQE2 BWF 
	BiologicalWeightingFunction LQE2BWF;
	LQE2BWF.SetWeightingFunction([](double const* params, double linealEnergy) {return ((params[0]*linealEnergy+params[1]*linealEnergy*linealEnergy)*std::exp(-params[2]*linealEnergy*linealEnergy));}, 3);

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
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.SetCellStudyParameters(H460FittingParams);
	BWFFitter.Initialize();
	BWFFitter.SetBetaParameterLowerConstraint(5,0);
	BWFFitter.SetBetaParameterLowerConstraint(4,0);
	BWFFitter.SetBetaParameterLowerConstraint(3,0);
	BWFFitter.SetBetaParameterLowerConstraint(2,0);
	BWFFitter.SetBetaParameterLowerConstraint(1,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto linearResults = BWFFitter.Fit();
	linearResults.PrintBasic();
	linearResults.PrintAIC(83);

	std::cout << std::endl << "Quadratic" << std::endl;
	QuadraticBWF.SetValues(std::vector<double> {0, 0.0139333, 0.127776});
	BWFFitter.SetAlphaWeightingFunction(QuadraticBWF);
	FifthBWF.SetValues(std::vector<double> {0, 1.89075e-07, 0, 0, 0, 0.130403});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(5,0);
	BWFFitter.SetBetaParameterLowerConstraint(4,0);
	BWFFitter.SetBetaParameterLowerConstraint(3,0);
	BWFFitter.SetBetaParameterLowerConstraint(2,0);
	BWFFitter.SetBetaParameterLowerConstraint(1,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto quadraticResults = BWFFitter.Fit();
	quadraticResults.PrintBasic();
	quadraticResults.PrintAIC(83);

	std::cout << std::endl << "Cubic" << std::endl;
	CubicBWF.SetValues(std::vector<double> {0, 0.000873803, 0, 0.151149});
	BWFFitter.SetAlphaWeightingFunction(CubicBWF);
	FifthBWF.SetValues(std::vector<double> {0, 5.37076e-08, 3.61649e-07, 0, 0, 0.134576});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(5,0);
	BWFFitter.SetBetaParameterLowerConstraint(4,0);
	BWFFitter.SetBetaParameterLowerConstraint(3,0);
	BWFFitter.SetBetaParameterLowerConstraint(2,0);
	BWFFitter.SetBetaParameterLowerConstraint(1,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto cubicResults = BWFFitter.Fit();
	cubicResults.PrintBasic();
	cubicResults.PrintAIC(83);

	std::cout << std::endl << "Fourth order polynomial" << std::endl;
	FourthBWF.SetValues(std::vector<double> {0, 1.44958e-05, 0, 0, 0.160277});
	BWFFitter.SetAlphaWeightingFunction(FourthBWF);
	FifthBWF.SetValues(std::vector<double> {0, 6.2094e-08, 0, 0, 0, 0.135782});
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(5,0);
	BWFFitter.SetBetaParameterLowerConstraint(4,0);
	BWFFitter.SetBetaParameterLowerConstraint(3,0);
	BWFFitter.SetBetaParameterLowerConstraint(2,0);
	BWFFitter.SetBetaParameterLowerConstraint(1,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto fourthResults = BWFFitter.Fit();
	fourthResults.PrintBasic();
	fourthResults.PrintAIC(83);

	// The result this converges to, is literally the same as the fourth order result, the first parameter goes to zero
	std::cout << std::endl << "Fifth order polynomial" << std::endl;
	FifthBWF.SetValues(std::vector<double> {0, 0, 1.44943e-05, 0, 0, 0.160252});
	BWFFitter.SetAlphaWeightingFunction(FifthBWF);
	FifthBWF.SetValues(std::vector<double> { 0, 6.20455e-08, 0, 0, 0, 0.135764}); //Fifth BWF likes getting zeros for beta, the others don't
	BWFFitter.SetBetaWeightingFunction(FifthBWF);
	BWFFitter.Initialize();
	BWFFitter.SetAlphaParameterLowerConstraint(5,0);
	BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	BWFFitter.SetBetaParameterLowerConstraint(5,0);
	BWFFitter.SetBetaParameterLowerConstraint(4,0);
	BWFFitter.SetBetaParameterLowerConstraint(3,0);
	BWFFitter.SetBetaParameterLowerConstraint(2,0);
	BWFFitter.SetBetaParameterLowerConstraint(1,0);
	BWFFitter.SetBetaParameterLowerConstraint(0,0);
	auto fifthResults = BWFFitter.Fit();
	fifthResults.PrintBasic();
	fifthResults.PrintAIC(83);

		//
		// Mairani et al.
		//

	// std::cout << std::endl << "Q" << std::endl;
	// QBWF.SetValues(std::vector<double> {0.0001});
	// BWFFitter.SetAlphaWeightingFunction(QBWF);
	// // FifthBWF.SetValues(std::vector<double> {0.1,0.1,0.1,0.1,0.1,0.1});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto QResults = BWFFitter.Fit();
	// QResults.PrintBasic();
	// QResults.PrintAIC(83);

	// std::cout << std::endl << "LE" << std::endl;
	// LEBWF.SetValues(std::vector<double> {0., 0.0139333});
	// BWFFitter.SetAlphaWeightingFunction(LEBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto LEResults = BWFFitter.Fit();
	// LEResults.PrintBasic();
	// LEResults.PrintAIC(83);

	// std::cout << std::endl << "QE" << std::endl;
	// QEBWF.SetValues(std::vector<double> {0.0, 0.000873803});
	// BWFFitter.SetAlphaWeightingFunction(QEBWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(1,-2);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,1e-6);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto QEResults = BWFFitter.Fit();
	// QEResults.PrintBasic();
	// QEResults.PrintAIC(83);

	// std::cout << std::endl << "LQE" << std::endl;
	// LQEBWF.SetValues(std::vector<double> {0.0, 0.000873803, 0.0139333});
	// BWFFitter.SetAlphaWeightingFunction(LQEBWF); //The ordering is strange, you have to specify the function first
	// BWFFitter.Initialize(); //then initialize
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0); //only then can you give the constraints
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto LQEResults = BWFFitter.Fit();
	// LQEResults.PrintBasic();
	// LQEResults.PrintAIC(83);

	// std::cout << std::endl << "LE2" << std::endl;
	// LE2BWF.SetValues(std::vector<double> {0., 0.0139333});
	// BWFFitter.SetAlphaWeightingFunction(LE2BWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto LE2Results = BWFFitter.Fit();
	// LE2Results.PrintBasic();
	// LE2Results.PrintAIC(83);

	// std::cout << std::endl << "QE2" << std::endl;
	// QE2BWF.SetValues(std::vector<double> {0.0, 0.000873803});
	// BWFFitter.SetAlphaWeightingFunction(QE2BWF);
	// BWFFitter.Initialize();
	// // BWFFitter.SetAlphaParameterLowerConstraint(1,-0.0005);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto QE2Results = BWFFitter.Fit();
	// QE2Results.PrintBasic();
	// QE2Results.PrintAIC(83);

	// std::cout << std::endl << "LQE2" << std::endl;
	// LQE2BWF.SetValues(std::vector<double> {0.0, 0.000873803, 0.0139333});
	// BWFFitter.SetAlphaWeightingFunction(LQE2BWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto LQE2Results = BWFFitter.Fit();
	// LQE2Results.PrintBasic();
	// LQE2Results.PrintAIC(83);

		//
		//Fitting the other functional forms
		//

	// std::cout << "Morstin" << std::endl;
	// MorstinBWF.SetValues(std::vector<double> {2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000 });
	// BWFFitter.SetAlphaWeightingFunction(MorstinBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// auto MorstinResults = BWFFitter.Fit();
	// MorstinResults.PrintBasic();
	// MorstinResults.PrintAIC(83);

	// std::cout << "Morstin plus offset" << std::endl;
	// MorstinPlusOffsetBWF.SetValues(std::vector<double> {0.0, 2.*std::pow(10,-7), 2.1*std::pow(10,-5), 2.5*std::pow(10,-6), 11460.000000});
	// BWFFitter.SetAlphaWeightingFunction(MorstinPlusOffsetBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(4,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// auto MorstinOffsetResults = BWFFitter.Fit();
	// MorstinOffsetResults.PrintBasic();
	// MorstinOffsetResults.PrintAIC(83);

	// std::cout << "Gaussian" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(2,1);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,1);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,1);
	// // BWFFitter.SetAlphaParameterUpperConstraint(1,500); 	//Overkill shouldn't happen any later than this.
	// // BWFFitter.SetAlphaParameterUpperConstraint(0,10000); 	//Prevent runaway amplitude
	// auto GaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianVarAmplitudeResults.PrintBasic();
	// GaussianVarAmplitudeResults.PrintAIC(83);

	// std::cout << "Gaussian w Gaussian Beta" << std::endl;
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {111.38, 517.436, 6975.15});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudeBWF);
	// GaussianVariableAmplitudeBWF.SetValues(std::vector<double> {415.535, 301.313, 0.184074});
	// BWFFitter.SetBetaWeightingFunction(GaussianVariableAmplitudeBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(1,500); 	//Overkill shouldn't happen any later than this.
	// BWFFitter.SetAlphaParameterUpperConstraint(0,10000); 	//Prevent runaway amplitude
	// auto GaussianGaussianVarAmplitudeResults = BWFFitter.Fit();
	// GaussianGaussianVarAmplitudeResults.PrintBasic();
	// GaussianGaussianVarAmplitudeResults.PrintAIC(83);

	// Performs just as bad as the old thing
	// std::cout << "Gaussian, new thing" << std::endl;
	// GaussianPlusOffsetBWF.SetValues(std::vector<double> {40, 80, 3});
	// BWFFitter.SetAlphaWeightingFunction(GaussianPlusOffsetBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(0,1000); 	//Prevent runaway amplitude
	// auto GaussianNewResults = BWFFitter.Fit();
	// GaussianNewResults.PrintBasic();

	//Not really getting good or interesting results
	// std::cout << "Gaussian plus offset" << std::endl;
	// GaussianVariableAmplitudePlusOffsetBWF.SetValues(std::vector<double> {1, 100, 100, 30});
	// BWFFitter.SetAlphaWeightingFunction(GaussianVariableAmplitudePlusOffsetBWF);
	// FifthBWF.SetValues(std::vector<double> {0, 0, 0, 0, 0, 0});
	// BWFFitter.SetBetaWeightingFunction(FifthBWF);
	// BWFFitter.Initialize();
	// //Precent runaway amplitude
	// BWFFitter.SetAlphaParameterLowerConstraint(3,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(2,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(1,0);
	// BWFFitter.SetAlphaParameterLowerConstraint(0,0);
	// BWFFitter.SetAlphaParameterUpperConstraint(0,1000);
	// // BWFFitter.SetBetaParameterLowerConstraint(5,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(4,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(3,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(2,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(1,0);
	// // BWFFitter.SetBetaParameterLowerConstraint(0,0);
	// auto GaussianVarAmplitudePlusOffsetResults = BWFFitter.Fit();
	// GaussianVarAmplitudePlusOffsetResults.PrintBasic();

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
	// BWFFunctionPlotter(c, legend, lineStyle, "", FifthBWF, fifthResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	// lineStyle.SetLineColor(kPink+2);
	BWFFunctionPlotter(c, legend, lineStyle, "", CubicBWF, cubicResults.alphaFunc.GetFittingParams(), "AL", 0., 100.);
	lineStyle.SetLineColor(kBlue+2);
	BWFFunctionPlotter(c, legend, lineStyle, "", FourthBWF, fourthResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	lineStyle.SetLineColor(kRed+2);
	BWFFunctionPlotter(c, legend, lineStyle, "", LinearBWF, linearResults.alphaFunc.GetFittingParams(), "L", 0., 100.);
	lineStyle.SetLineColor(kOrange+2);
	BWFFunctionPlotter(c, legend, lineStyle, "", QuadraticBWF, quadraticResults.alphaFunc.GetFittingParams(), "L", 0., 100.);

	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/polynomialhalpha_CERES_BWFs_betaFifth.jpg";
	c->SaveAs((TString)outputName); 

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

void Survival_Fitting()
{
	SetupH460SurvivalParameters();

	AlphaBeta_Fitter fitter{};
	fitter.SetCellStudyParameters(H460FittingParams);
	double* AlphaBeta = fitter.Fit(nullptr,false);

	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.95,0.72+0.23);//x start, y start, x end, yend
	// auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
	legend->SetTextSize(0.05);

	TAttMarker markerAtts;
	markerAtts.SetMarkerColor(kRed+2);
	markerAtts.SetMarkerSize(8);
	markerAtts.SetMarkerStyle(8);

	PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, false);
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/alpha_H460.jpg";
	c->SaveAs((TString)outputName); 

	PlotAlphaBeta(c, legend, "", markerAtts, "AP", H460FittingParams, AlphaBeta, true);
	outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/beta_H460.jpg";
	c->SaveAs((TString)outputName); 
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

void H460_Ceres()
{
	Fifth_Beta();
	// Survival_Fitting();
}