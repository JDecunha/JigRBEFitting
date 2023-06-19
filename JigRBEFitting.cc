//Fitting specific cell types
#include "H460Fitting.h"
#include "LinealSpectra.h"
#include "ProtonSpectra.h"
#include "Utilities.h"
#include "LETCalculator.h"
#include "H460_Ceres.h"
#include "H460_LET_Ceres.h"
#include "H1437_Ceres.h"
#include "H1437_LET_Ceres.h"
#include "Uwes_Spectrum.h"
#include "Spectra_For_Erik.h"

int main()
{
	// LETCalcMain();
	// Uwes_Spectrum();
	// H1437_LET_Ceres();
	// LinealSpectra::SaveKeWeightedLinealSpectra();
	// auto lib = LinealSpectra::GetKeWeightedLinealSpectraJuly2023("1e3");
	// LinealSpectra::SaveKeWeightedLinealSpectraCSV();

	// auto lib = ErikSpectra::GetMonoenergeticLinealSpectra();
	// ErikSpectra::PlotMonoenergeticLinealSpectra(lib);
	// ErikSpectra::PrintLESSpectra(lib);

	H460_Ceres();
	
};