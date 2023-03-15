#include "BiologicalWeightingFunction.h"

void BiologicalWeightingFunction::SetWeightingFunction(double (*BWF)(double const*, double), int nParams)
{
	_pBWF = BWF;
	_nFittingParams = nParams;
	_params = new double [_nFittingParams](); //parentheses are important, there zero initialize the array
}
