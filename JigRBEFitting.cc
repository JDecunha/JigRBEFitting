//Fitting specific cell types
#include "H460Fitting.h"
#include "LinealSpectra.h"
#include "Utilities.h"
#include "LETCalculator.h"
#include "H460_Ceres.h"

int main()
{
	H460_Ceres();
	H460Fitting();
	// LETCalcMain();

	// gInterpreter->GenerateDictionary("pair<string,TH1D>;vector<pair<string,TH1D> >", "TH1.h;string;utility;vector");
	// auto keWeightedSpectra = LinealSpectra::GetKeWeightedLinealSpectra("1e3");
	// LinealSpectra::PlotKeWeightedLinealSpectraMultigraph(keWeightedSpectra);
};