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
#include "ProtonSpectra.h" //Proton KESpectra
#include "LinealSpectra.h" //Lineal energy spectra from KE Spectra
//Functions for fitting and plotting
#include "SurvivalPlotting.h"
#include "BWF_Fitter.h"
#include "BWF_Fitter_Beta.h"
#include "BWF_Plotter.h"
//This file
#include "H460Fitting.h"

//I put the survival parameters at global scope
CellStudySurvivalParameters H460Params{};

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

	//Photon beta (from Fada's calculations)
	double beta = 0.097;

	//Manually fitted alpha-beta for each experiment (from Fada's calculations)
	std::vector<double> alphas {0.268, 0.226, 0.151, 0.150, 0.166, 0.137, 0.206, 0.117, 0.318, 0.446, 0.596, 0.883};
	std::vector<double> betas {0.097, 0.112, 0.134, 0.134, 0.134, 0.146, 0.125, 0.159, 0.154, 0.341, 0.662, 0.956};

	//Get the lineal energy spectra library
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");
	std::vector<std::pair<std::string, TH1D>>* dySpectraLibrary;
	TFile* dySpectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraCellStudy.root");
	dySpectrumFile->GetObject("Lineal_energy_library",dySpectraLibrary);
	dySpectrumFile->Close(); //Close file

	//Set the cell study parameters
	H460Params.dySpectra = *dySpectraLibrary;
	H460Params.dose = doses;
	H460Params.survivingFraction = survivingFractions;
	H460Params.survivingFractionUncertainty = survivingFractionUncertainty;
	H460Params.beta = beta;	
	H460Params.LETd = LETd;
}

	
void H460Fitting() 
{
	//
	//Fitting
	//

	//This function takes the H460 cell study params that are at global scope and fills them with the dose, SF, and d(y) information.
	SetupH460SurvivalParameters();

	//Set up the fitter
	BWF_Fitter_Beta fitter{};
	fitter.SetCellStudyParameters(H460Params);
	BWF_Fitter fitterold{};
	fitterold.SetCellStudyParameters(H460Params);


	//Create linear BWF
	BiologicalWeightingFunction LinearandFixedBWF;
	LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);
	//Linear function fitting
	fitter.SetWeightingFunction(LinearandFixedBWF);
	double linearInitialGuess [] = {1, 1, 0.097};
	double* LinearandFixedParams = fitter.Fit(linearInitialGuess,true);

	//Create quadratic BWF
	BiologicalWeightingFunction QuadraticLinearFixedBWF;
	QuadraticLinearFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[2]*linealEnergy*linealEnergy+params[1]*linealEnergy+params[0]);}, 3);
	//BWF Function Fitting
	fitterold.SetWeightingFunction(QuadraticLinearFixedBWF);
	double quadraticInitialGuess [] = {0.3, 0.1, 0};//, 0.097};
	double* QuadraticLinearFixedParams = fitterold.Fit(quadraticInitialGuess,true);
	fitter.SetWeightingFunction(QuadraticLinearFixedBWF);
	double quadraticInitialGuess2 [] = {0.3, 0.1, 0, 0.097};//, 0.097};
	double* QuadraticLinearFixedParams2 = fitter.Fit(quadraticInitialGuess2,true);


	//
	//Plotting
	//
	
	//Plotting the survival data

	//Setup the canvas
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	auto legend = new TLegend(0.52,0.72,0.89,0.72+0.16);//x start, y start, x end, yend
	legend->SetTextSize(0.05);

	//Setup the marker attributes
	TAttLine lineStyle{};
	lineStyle.SetLineColor(kGreen+2);
	lineStyle.SetLineWidth(3);
	lineStyle.SetLineStyle(2);

	//Plot the survival data
	SurvivalDataMultigraph(c, legend, H460Params);

	//Plot the BWF prediction on the survival data
	lineStyle.SetLineColor(kRed+2);
	GeneralizedBWFMultigraphPlotterBeta(c, legend, lineStyle, "Linear BWF", H460Params, LinearandFixedBWF, LinearandFixedParams);

	lineStyle.SetLineColor(kBlue+4);
	GeneralizedBWFMultigraphPlotterBeta(c, legend, lineStyle, "Quadratic BWF", H460Params, QuadraticLinearFixedBWF, QuadraticLinearFixedParams2);
	


	//Save
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/FitswithBeta.jpg";
	c->SaveAs((TString)outputName); 
























	//Set up the fitter
	// BWF_Fitter fitterold{};
	// fitterold.SetCellStudyParameters(H460Params);

	// //Linear function fitting
	// BiologicalWeightingFunction LinearandFixedBWF;
	//LinearandFixedBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (params[1]*linealEnergy)+params[0];}, 2);
	// fitterold.SetWeightingFunction(LinearandFixedBWF);
	// double* LinearandFixedParams2 = fitterold.Fit(nullptr,true);



	// //Gaussian function fitting
	// BiologicalWeightingFunction GaussianBWF;
	// GaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
	// { 
	// 	return ((6000/(params[0]))*std::exp(-((linealEnergy-params[1])*(linealEnergy-params[1]))/(params[0]*params[0]*2))); 
	// }
	// , 2);
	// fitter.SetWeightingFunction(GaussianBWF);
	// double gaussianInitialGuess [] = {70, 3};
	// double* GaussianParams = fitter.Fit(gaussianInitialGuess,true);

	// //Skew function fitting
	// BiologicalWeightingFunction SkewGaussianBWF;
	// SkewGaussianBWF.SetWeightingFunction( [] (double* params, double linealEnergy) 
	// { 
	// 	//With regards to this wonderful man on the ROOT forums: https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	//     double xi = params[0];
	//     double omega = params[1];
	//     double alpha = params[2];
	//     double arg = (linealEnergy - xi) / omega;
	//     double smallphi = TMath::Gaus(arg, 0.0, 1.0, true);
	//     double bigphi = 0.5 * (1 + std::erf(alpha * arg/std::sqrt(2)));
	//     return (2200)/omega * smallphi * bigphi;
	// }
	// , 3);
	// fitter.SetWeightingFunction(SkewGaussianBWF);
	// double skewGaussianInitialGuess [] = {115, 40, -8};
	// double* skewGaussianParams = fitter.Fit(skewGaussianInitialGuess,true);


	//double skewGaussianParams[4];
	//skewGaussianParams[0] = 204; skewGaussianParams[1] = 58; skewGaussianParams[2] = 0; skewGaussianParams[3] = 14000;

	//On creating a skew gaussian in C++
	//https://root-forum.cern.ch/t/how-to-fit-a-skew-gaussian/50922/2
	//

	//Morstin (1989).
	/*BiologicalWeightingFunction MorstinBWF;
	MorstinBWF.SetWeightingFunction([](double* params, double linealEnergy) {return (1-std::exp(linealEnergy*params[0]+params[1]*linealEnergy*linealEnergy));}, 2);
	fitter.SetWeightingFunction(MorstinBWF);
	double morstinInitialGuess [] = {0.1, 5};
	fitter.Fit(morstinInitialGuess);*/
	

	

	// lineStyle.SetLineColor(kRed+2);

	// GeneralizedBWFMultigraphPlotter(c, legend, lineStyle, "Gaussian BWF", H460Params, GaussianBWF, GaussianParams);

	// lineStyle.SetLineColor(kRed+2);
	// GeneralizedBWFMultigraphPlotter(c, legend, lineStyle, "Skew-Gaussian BWF", H460Params, SkewGaussianBWF, skewGaussianParams);
	//GeneralizedBWFMultigraphPlotter(c, legend, lineStyle, "Mortsin BWF", H460Params, MorstinBWF, MorstinParams);



	//Plotting the BWFs themselves

	//Setup the canvas
	// gStyle->SetOptStat(0); //Don't print the stats window in the top right
	// TCanvas* c = new TCanvas("c","c");
	// c->SetCanvasSize(4000, 2400);

	// auto legend = new TLegend(0.14,0.72,0.37,0.72+0.16);//x start, y start, x end, yend
	// legend->SetTextSize(0.04);

	// //Setup the marker attributes
	// TAttLine lineStyle{};
	
	// lineStyle.SetLineWidth(5);
	// lineStyle.SetLineStyle(1);

	// //Plot the BWFs
	// lineStyle.SetLineColor(kBlue+4);
	// BWFFunctionPlotter(c, legend, lineStyle, "Gaussian BWF", GaussianBWF, GaussianParams,"AL");

	// lineStyle.SetLineColor(kGreen+3);
	// //BWFFunctionPlotter(c, legend, lineStyle, "Linear BWF", LinearandFixedBWF, LinearandFixedParams ,"L");

	// lineStyle.SetLineColor(kRed+2);
	// //BWFFunctionPlotter(c, legend, lineStyle, "Quadratic BWF",QuadraticLinearFixedBWF, QuadraticLinearFixedParams,"L");

	// lineStyle.SetLineColor(kAzure+2);
	// BWFFunctionPlotter(c, legend, lineStyle, "Skew Gaussian BWF",SkewGaussianBWF, skewGaussianParams,"L");

	// lineStyle.SetLineColor(kRed);
	// BWFFunctionPlotter(c, legend, lineStyle, "Skew Gaussian Initial Guess",SkewGaussianBWF, skewGaussianInitialGuess,"L");

	// //GeneralizedBWFMultigraphPlotter(c, legend, lineStyle, "Mortsin BWF", H460Params, MorstinBWF, MorstinParams);

	// //Save
	// std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/BWFsNovember3rd.jpg";
	// c->SaveAs((TString)outputName); 

}
