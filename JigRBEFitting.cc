//Fitting specific cell types
#include "H460Fitting.h"
#include "LinealSpectra.h"
#include "Utilities.h"
#include "LETCalculator.h"
#include "H460_Ceres.h"
#include "H460_LET_Ceres.h"
#include "H1437_Ceres.h"
#include "H1437_LET_Ceres.h"
#include "Uwes_Spectrum.h"

int main()
{
	// auto out = LinealSpectra::GetKeWeightedFrequencyLinealSpectra("1e3");
	// LinealSpectra::PlotKeWeightedLinealSpectra(out);
	// LinealSpectra::SaveKeWeightedFrequencyLinealSpectra();
	Uwes_Spectrum();
};