#include "BiologicalWeightingFunction.h"

void BiologicalWeightingFunction::SetWeightingFunction(double (*BWF)(double*, double), int nParams)
{
	_pBWF = BWF;
	_nFittingParams = nParams;
}
