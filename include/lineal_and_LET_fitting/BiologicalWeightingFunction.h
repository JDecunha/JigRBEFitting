#pragma once

#include <algorithm>
#include <iostream>
#include <vector>



class BiologicalWeightingFunction 
{
	public:	
		void SetWeightingFunction(double (*BWF)(double const*, double), int nParams);
		void SetValues(std::vector<double> values) 
		{ 
			if(values.size() == _nFittingParams)  { std::reverse(values.begin(),values.end()); _params = values; } //We call reverse most signifigant fitting param is first
			else {std::cout << "Error: Number of values does not match number of fitting params." << std::endl; throw;}
		};

		double GetValue(double const* inputs, double linealEnergy) const { return _pBWF(inputs,linealEnergy); };
		double GetValue(double linealEnergy) const { return _pBWF(_params.data(),linealEnergy); };
		double* GetFittingParams() { return _params.data(); }
		int GetNumFittingParams() const { return _nFittingParams; };

		void PrintFitParams() 
		{
			for (int i = _params.size()-1; i >= 0; --i) 
			{
				std::cout << _params[i];
				if (i != 0) {std::cout << ", ";}
			} 
		}

	private:
		double (*_pBWF)(double const*, double);
		int _nFittingParams;
		std::vector<double> _params;
};