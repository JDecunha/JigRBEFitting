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

	//Plot
	//Setup the canvas
	gStyle->SetOptStat(0); //Don't print the stats window in the top right
	TCanvas* c = new TCanvas("c","c");
	c->SetCanvasSize(9000, 5000);
	c->SetFillStyle(4000);
	c->SetFrameFillStyle(4000);
	c->Divide(4,3,0.000000005,0.001);
	//auto legend = new TLegend(0.52,0.72,0.89,0.72+0.16);//x start, y start, x end, yend
	auto legend = new TLegend(0.14,0.14,0.14+0.33,0.14+0.16);
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
	GeneralizedBWFMultigraphPlotterBeta(c, legend, lineStyle, "Linear, 1 um", H460Params, LinearandFixedBWF, LinearandFixedParams);

	//Get the lineal energy spectra library
	gInterpreter->GenerateDictionary("pair<string,TH1F>;vector<pair<string,TH1F> >", "TH1.h;string;utility;vector");
	std::vector<std::pair<std::string, TH1D>>* dySpectraLibrary;
	TFile* dySpectrumFile = TFile::Open("/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Data/LinealSpectraCellStudy_100nm.root");
	dySpectrumFile->GetObject("Lineal_energy_library",dySpectraLibrary);
	dySpectrumFile->Close(); //Close file

	// //Re-fit with 100 nm microdosimetric spectra
	H460Params.dySpectra = *dySpectraLibrary;
	fitter.SetCellStudyParameters(H460Params);
	double linearInitialGuess_100nm [] = {0, 0.03, 0.097};
	double* LinearandFixedParams_100nm = fitter.Fit(linearInitialGuess_100nm,true);

	//Plotting the survival data
	lineStyle.SetLineColor(kBlue+4);
	GeneralizedBWFMultigraphPlotterBeta(c, legend, lineStyle, "Linear, 100 nm", H460Params, LinearandFixedBWF, LinearandFixedParams_100nm);
	
	//Save
	std::string outputName = "/home/joseph/Dropbox/Documents/Work/Projects/MDA_vitro_RBE/Images/fitting/100nm_fit_test.jpg";
	c->SaveAs((TString)outputName); 