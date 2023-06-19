#pragma once

std::pair<std::vector<double>,std::vector<double>> PullStoppingPowersFromCsv();
std::vector<std::pair<float,float>> ComputeUwesLETs();
void ComputeLETs();
void LETCalcMain();