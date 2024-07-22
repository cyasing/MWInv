#pragma once
#include "OptimizeBase.h"

// Newton Solver Based optimization class using dlib
class OptimizeNW :
	public OptimizeBase
{
public:
	OptimizeNW();
	~OptimizeNW();
	void solve_optimization_problem();
};

