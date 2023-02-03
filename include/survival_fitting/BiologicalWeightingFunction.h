#pragma once

class BiologicalWeightingFunction 
{
	public:	
		void SetWeightingFunction(double (*BWF)(double*, double), int nParams);
		double GetValue(double* inputs, double linealEnergy) const { return _pBWF(inputs,linealEnergy); };
		int GetNumFittingParams() const { return _nFittingParams; };

	private:
		double (*_pBWF)(double*, double);
		int _nFittingParams;
};