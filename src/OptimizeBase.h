#pragma once
#include "OptimizationParams.h"

//Base class for optimization algorithms
class OptimizeBase
{
private:
	bool print_inv;
	bool save_inv;
public:
	OptimizeBase();
	virtual ~OptimizeBase();
	virtual void solve_optimization_problem(bool print_inv, bool save_inv) = 0;
};