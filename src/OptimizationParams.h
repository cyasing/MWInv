#pragma once
#include <string>
#include "dlib/optimization.h" //may need more headers

// Parameters for the optimization methods
// Made this struct since different methods may require different parameters
// Will use certain parameters from this for certain methods
struct OptimizationParams
{
public:
	std::string method = "cg";
	float min_error = 1e-5;
};

OptimizationParams::OptimizationParams(const OptimizationParams& opt_params)
{
	this->method = opt_params.method;
	this->min_error = opt_params.min_error;
}