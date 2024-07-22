#pragma once
#include "OptimizeBase.h"

// Conjugate Gradient optimization class using dlib
class OptimizeCG :
    public OptimizeBase
{
public:
	OptimizeCG();
	~OptimizeCG();
	void solve_optimization_problem();
};

