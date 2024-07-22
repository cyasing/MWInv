#pragma once

#include "MWInvSystemSetup.h"
#include "DefaultMesh.h"
#include "MFEM/mfem.hpp"
#include <functional>


// MFEM Custom Type Declaration
typedef double real_t;
template <typename T> T pow2(const T& x) { return x * x; }

// Global Variables for the Forward Solver
mfem::Array2D<real_t> comp_domain_bdr;
mfem::Array2D<real_t> domain_bdr;

// Global Variables for the Forward Solver
real_t mu = 1.0;
real_t epsilon = 1.0;
real_t omega;
int dim;

// Class for calling the forward solver of MFEM based on system and mesh setup
class ForwardSolve
{
private:
	// System and Mesh Objects
	MWInvSystemSetup& system;
	DefaultMesh& region_mesh;
private:
	/* MFEM Variables */

	// Args Parser: Default Inputs  
    const char* mesh_file = nullptr;
    int order = 1;
    int ref_levels = 3;
    real_t freq = 4.0;
    bool herm_conv = true;
    bool umf_solver = false;
    bool visualization = 1;
    bool pa = false;
    const char* device_config = "cpu";

	// Region Variables
    mfem::Array2D<real_t> length;
    mfem::Mesh* mesh;
    mfem::Array2D<real_t> comp_domain_bdr;
    mfem::Array<int> ess_tdof_list;
    mfem::Array<int> ess_bdr;
    mfem::Array<int> attr;
    mfem::Array<int> attrPML;

	// Physical Constants
    mfem::ConstantCoefficient muinv;
    mfem::ConstantCoefficient omeg;


	// FE Variables
    mfem::FiniteElementCollection* fec;
    mfem::FiniteElementSpace* fespace;
    mfem::FiniteElementSpace* fespace_x;
    PML* pml;
    mfem::ComplexOperator::Convention conv;

    mfem::VectorFunctionCoefficient f;
    mfem::ComplexGridFunction x;
    mfem::ComplexLinearForm b;
    mfem::SesquilinearForm a;

    std::vector<mfem::Vector> e_ref_real_rx;
    std::vector<mfem::Vector> e_ref_imag_rx;

    mfem::Vector circle_center;
    mfem::FunctionCoefficient contrast_coeff;
    IterContrastCoefficient sq_omega_eps;

    mfem::RestrictedCoefficient restr_muinv;
    mfem::RestrictedCoefficient restr_omeg;
    mfem::RestrictedCoefficient restr_sq_omega_eps;
    mfem::FunctionCoefficient contrast_coeff;

public:
	// FEM Results for the inverse problem
	mfem::Vector b_inv;
	mfem::SparseMatrix A_inv;
public:
    ForwardSolve(MWInvSystemSetup& system, DefaultMesh& mesh);
	~ForwardSolve();
    int forward_solve_setup(int argc, char* argv[]);
	void forward_solve_run();
	void set_inverse_system();
	double evaluate_residual();
	mfem::Vector evaluate_residual_derivative();

private:

    bool IsWithinSquare(const mfem::Vector& point, const mfem::Vector& center, double half_size);
    
    void eval_source_receivers_locations(
        const size_t num_antennas,
        const double radius,
        const size_t iteration,
        Points2D<double>& location_source,
        Points2D<double>& location_receivers);

	void source(const mfem::Vector& x, mfem::Vector& f);
};

// Iterating Contrast Coefficient
class IterContrastCoefficient : public mfem::Coefficient {
public:
    IterContrastCoefficient(
        real_t omega_,
        real_t default_epsilon,
        real_t anomaly_epsilon,
        const mfem::Vector& center,
        real_t radius)
        :
        omega(omega_),
        default_epsilon(default_epsilon),
        circle_epsilon(anomaly_epsilon),
        circle_center(center),
        circle_radius(radius) {}

    // Necessary when deriving from Coefficient
    virtual double Eval(mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override
    {
        mfem::Vector x;
        T.Transform(ip, x);
        return Evaluate_Coefficient(x);
    }

    double Evaluate_Coefficient(const mfem::Vector& x)
    {
        if (pow2(x[0] - circle_center[0]) + pow2(x[1] - circle_center[1]) <= pow2(circle_radius))
        {
            return circle_epsilon * pow2(omega);
        }
        else
        {
            return default_epsilon * pow2(omega);
        }
    }
private:
    real_t omega;
    real_t default_epsilon;
    real_t circle_epsilon;
    mfem::Vector circle_center;
    real_t circle_radius;
};

// Class for setting up a simple Cartesian PML region
class PML
{
private:
    mfem::Mesh* mesh;

    int dim;

    // Length of the PML Region in each direction
    mfem::Array2D<real_t> length;

    // Computational Domain Boundary
    mfem::Array2D<real_t> comp_dom_bdr;

    // Domain Boundary
    mfem::Array2D<real_t> dom_bdr;

    // Integer Array identifying elements in the PML
    // 0: in the PML, 1: not in the PML
    mfem::Array<int> elems;

    // Compute Domain and Computational Domain Boundaries
    void SetBoundaries();

public:
    // Constructor
    PML(mfem::Mesh* mesh_, mfem::Array2D<real_t> length_);

    // Return Computational Domain Boundary
    mfem::Array2D<real_t> GetCompDomainBdr() { return comp_dom_bdr; }

    // Return Domain Boundary
    mfem::Array2D<real_t> GetDomainBdr() { return dom_bdr; }

    // Return Markers list for elements
    mfem::Array<int>* GetMarkedPMLElements() { return &elems; }

    // Mark elements in the PML region
    void SetAttributes(mfem::Mesh* mesh_);

    // PML complex stretching function
    void StretchFunction(const mfem::Vector& x, std::vector<std::complex<real_t>>& dxs);
};

// Class for returning the PML coefficients of the bilinear form
class PMLDiagMatrixCoefficient : public mfem::VectorCoefficient
{
private:
    PML* pml = nullptr;
    void (*Function)(const mfem::Vector&, PML*, mfem::Vector&);
public:
    PMLDiagMatrixCoefficient(int dim, void(*F)(const mfem::Vector&, PML*,
        mfem::Vector&),
        PML* pml_)
        : VectorCoefficient(dim), pml(pml_), Function(F)
    {}

    using VectorCoefficient::Eval;

    virtual void Eval(mfem::Vector& K, mfem::ElementTransformation& T,
        const mfem::IntegrationPoint& ip)
    {
        real_t x[3];
        mfem::Vector transip(x, 3);
        T.Transform(ip, transip);
        K.SetSize(vdim);
        (*Function)(transip, pml, K);
    }
};