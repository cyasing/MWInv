#pragma once
#include "ForwardSolve.h"
#include "OptimizeBase.h"

// Main class for the reconstruction of the model
// Unsure of whether to take inheritance or composition approach
class InverseSolve : public OptimizeBase, private ForwardSolve
{
private:
	OptimizationParams opt_params;
	static MWInvSystemSetup system;
	static DefaultMesh mesh;
	ForwardSolve em_scattering;
	IterContrastCoefficient contrast_coefficient;
public:
	InverseSolve();
	~InverseSolve();

	double forward_residual(const dlib::matrix<double>& m);
	dlib::matrix<double> residual_derivatives(const dlib::matrix<double>& m);

	void solve_optimization_problem(bool print_inv, bool save_inv);
};