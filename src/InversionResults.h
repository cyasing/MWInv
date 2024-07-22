#pragma once
#include <vector>
#include <complex>

// Results of the inversion process: Saved from InverseSolve
// Could add more parameters as needed
struct InversionResults
{
	std::vector<std::vector<std::complex<double>>> reconstructed_permittivity;
	std::vector<double> convergence_rate;
	std::vector<double> error_rate;
};