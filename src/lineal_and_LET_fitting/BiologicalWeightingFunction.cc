#include "BiologicalWeightingFunction.h"

void BiologicalWeightingFunction::SetWeightingFunction(double (*BWF)(double const*, double), int nParams)
{
	_pBWF = BWF;
	_nFittingParams = nParams;

	for (int i = 0; i < nParams; ++i)
	{
		_params.push_back(0);
	}
}
