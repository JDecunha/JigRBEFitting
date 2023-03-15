#pragma once

class BiologicalWeightingFunction 
{
	public:	
		void SetWeightingFunction(double (*BWF)(double const*, double), int nParams);
		double GetValue(double const* inputs, double linealEnergy) const { return _pBWF(inputs,linealEnergy); };
		double *GetFittingParams() { return _params; }
		int GetNumFittingParams() const { return _nFittingParams; };

	private:
		double (*_pBWF)(double const*, double);
		int _nFittingParams;
		double* _params;
};