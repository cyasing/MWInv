#pragma once
#include "OptimizeBase.h"

// Least Squares optimization class using dlib
class OptimizeLS :
	public OptimizeBase
{
public:
	OptimizeLS();
	~OptimizeLS();
	void solve_optimization_problem();
};

