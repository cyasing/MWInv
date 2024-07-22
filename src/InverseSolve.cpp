#include "InverseSolve.h"

InverseSolve::InverseSolve()
	: em_scattering(system, mesh), opt_params()
{}

InverseSolve::~InverseSolve()
{}

double InverseSolve::forward_residual(const dlib::matrix<double>& m) // Input either GridFunction or MatrixCoefficient or just permittivity?
{
	// Set m to the permittivity guess value to the forward class IterContrastCoefficient
	mfem::Vector centre(2);
	centre[0] = 0.0;
	centre[1] = 0.0;
	contrast_coefficient = IterContrastCoefficient(omega, epsilon, m(0), centre, 1.0);
	// Set the permittivity values from input m using constructor
	em_scattering.forward_solve_run();
	em_scattering.set_inverse_system();
    double residual = em_scattering.evaluate_residual();
}

dlib::matrix<double> InverseSolve::residual_derivatives(const dlib::matrix<double>& m)
{
	mfem::Vector derivative = em_scattering.evaluate_residual_derivative();
	dlib::matrix<double> dlib_derivative(derivative.Size(), 1);
	for (int i = 0; i < derivative.Size(); i++)
	{
		dlib_derivative(i, 0) = derivative[i];
	}
	return dlib_derivative;
}

void InverseSolve::solve_optimization_problem(bool print_inv, bool save_inv)
{
	// Should use the Optimization classes to solve this, but for now, just a placeholder
    dlib::matrix<double> starting_point = { 4, 8 };

    std::cout << "Running the Inverse Solver: Minimising the Residual" << std::endl;

    if (opt_params.method == "cg")
    {
        dlib::find_min(dlib::cg_search_strategy(),  // Use CG search algorithm
            dlib::objective_delta_stop_strategy(opt_params.min_error), // Stop when the change in rosen() is less than 1e-7
            forward_residual, residual_derivatives, starting_point, -1);

        std::cout << "CG Optimization, Minimum Residual Value:\n" << starting_point << std::endl;
    }
}