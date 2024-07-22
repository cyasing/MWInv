#include "ForwardSolve.h"

constexpr real_t operator""_r(long double v)
{
    return static_cast<real_t>(v);
}

constexpr real_t operator""_r(unsigned long long v)
{
    return static_cast<real_t>(v);
}

ForwardSolve::ForwardSolve(MWInvSystemSetup& refsystem, DefaultMesh& refmesh)
	: system(refsystem), region_mesh(refmesh), 
    f(dim, source), x(fespace), b(fespace,conv), a(fespace, conv),
	restr_muinv(muinv, attr), restr_omeg(omeg, attr), restr_sq_omega_eps(contrast_coeff, attr)
{
    mesh_file = "../data/inline-quad.mesh";
    mesh = new mfem::Mesh(mesh_file, 1, 1);
	dim = mesh->Dimension();
	length.SetSize(dim, 2); length = 0.25;
	fec = new mfem::ND_FECollection(order, dim);
	fespace = new mfem::FiniteElementSpace(mesh, fec);
    dim = mesh->Dimension();

    pml = new PML(mesh, length);
    comp_domain_bdr = pml->GetCompDomainBdr();

    fespace_x = x.FESpace();

	freq = 4.0_r;
	omega = 2.0_r * M_PI * freq;
	muinv = mfem::ConstantCoefficient(1.0_r / mu);
	omeg = mfem::ConstantCoefficient(- pow2(omega) * epsilon);

    circle_center = mfem::Vector(dim);
    circle_center(0) = 0.0;
    circle_center(1) = 0.0;
    IterContrastCoefficient sq_omega_eps(omega, epsilon, 2.0, circle_center, 0.3);
    mfem::FunctionCoefficient contrast_coeff([&sq_omega_eps](const mfem::Vector& x) { return sq_omega_eps.Evaluate_Coefficient(x); });

    // 12c. Restricting the coefficients based on whether they are in the PML region or not
    contrast_coeff = mfem::FunctionCoefficient([&sq_omega_eps](const mfem::Vector& x) { return sq_omega_eps.Evaluate_Coefficient(x); });
    restr_sq_omega_eps = mfem::RestrictedCoefficient(contrast_coeff, attr);

    // 12c. Restricting the coefficients based on whether they are in the PML region or not
    mfem::RestrictedCoefficient restr_muinv(muinv, attr);
    mfem::RestrictedCoefficient restr_omeg(omeg, attr);
    mfem::RestrictedCoefficient restr_time_derivative_term(contrast_coeff, attr);


}
ForwardSolve::~ForwardSolve()
{}

int ForwardSolve::forward_solve_setup(int argc, char* argv[])
{
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order",
        "Finite element order (polynomial degree).");
    args.AddOption(&ref_levels, "-ref", "--refinements",
        "Number of refinements");
    args.AddOption(&mu, "-mu", "--permeability",
        "Permeability of free space (or 1/(spring constant)).");
    args.AddOption(&epsilon, "-eps", "--permittivity",
        "Permittivity of free space (or mass constant).");
    args.AddOption(&freq, "-f", "--frequency",
        "Frequency (in Hz).");
    args.AddOption(&herm_conv, "-herm", "--hermitian", "-no-herm",
        "--no-hermitian", "Use convention for Hermitian operators.");

    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
        "--no-visualization",
        "Enable or disable GLVis visualization.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
        "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&device_config, "-d", "--device",
        "Device configuration string, see Device::Configure().");
    args.Parse();

    // Enable hardware devices such as GPUs, and programming models such as
    //    CUDA, OCCA, RAJA and OpenMP based on command line options.
    mfem::Device device(device_config);
    device.Print();

    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);

    omega = real_t(2.0 * M_PI) * freq;

    // Refine the mesh to increase the resolution.
    for (int l = 0; l < ref_levels; l++)
    {
        mesh->UniformRefinement();
    }

    // Set element attributes in order to distinguish elements in the
    //    PML region
    pml->SetAttributes(mesh);

    // Define the Nedelec basis finite element space on the mesh.
    int size = fespace->GetTrueVSize();

    std::cout << "Number of finite element unknowns: " << size << std::endl;

    // Determine the list of true essential boundary dofs. In this example,
    //    the boundary conditions are defined based on the specific mesh and the
    //    problem type.

    if (mesh->bdr_attributes.Size())
    {
        ess_bdr.SetSize(mesh->bdr_attributes.Max());
        ess_bdr = 1;
    }
    fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    conv = herm_conv ? mfem::ComplexOperator::HERMITIAN : mfem::ComplexOperator::BLOCK_SYMMETRIC;

	// SOURCE RECEIVER LOCATIONS AND FUNCTION SETUP - can be skipped when integrated with MWInvSystemSetup class

    real_t n = 5_r * omega * sqrt(epsilon * mu) / real_t(M_PI);
    real_t intensity = pow2(n) / real_t(M_PI);
    size_t iteration = 0;
    eval_source_receivers_locations(
        system.num_antennas, 
        system.radius, 
        system.iteration, 
        system.location_source, 
        system.location_receivers);
    return 0;
}

void ForwardSolve::forward_solve_run()
{
    b.AddDomainIntegrator(NULL, new mfem::VectorFEDomainLFIntegrator(f));
    b.Vector::operator=(0.0);
    b.Assemble();

    // Define the solution vector x as a complex finite element grid function.
    x = 0.0;
    mfem::VectorFunctionCoefficient E_Re(dim, E_bdr_data_Re);
    mfem::VectorFunctionCoefficient E_Im(dim, E_bdr_data_Im);
    x.ProjectBdrCoefficientTangent(E_Re, E_Im, ess_bdr);

    if (mesh->attributes.Size())
    {
        attr.SetSize(mesh->attributes.Max());
        attrPML.SetSize(mesh->attributes.Max());
        attr = 0;   attr[0] = 1;
        attrPML = 0;
        if (mesh->attributes.Max() > 1)
        {
            attrPML[1] = 1;
        }
    }

    // Will need to extend the operators to handle complex coefficients based on example

    mfem::CurlCurlIntegrator* k = new mfem::CurlCurlIntegrator(restr_muinv);
    a.AddDomainIntegrator(k, NULL);

    mfem::VectorFEMassIntegrator* m = new mfem::VectorFEMassIntegrator(restr_sq_omega_eps);
    a.AddDomainIntegrator(m, NULL);

    // 12e. Defining the PML coefficients
    int cdim = (dim == 2) ? 1 : dim;
    PMLDiagMatrixCoefficient pml_c1_Re(cdim, detJ_inv_JT_J_Re, pml);
    PMLDiagMatrixCoefficient pml_c1_Im(cdim, detJ_inv_JT_J_Im, pml);
    mfem::ScalarVectorProductCoefficient c1_Re(muinv, pml_c1_Re);
    mfem::ScalarVectorProductCoefficient c1_Im(muinv, pml_c1_Im);
    mfem::VectorRestrictedCoefficient restr_c1_Re(c1_Re, attrPML);
    mfem::VectorRestrictedCoefficient restr_c1_Im(c1_Im, attrPML);

    PMLDiagMatrixCoefficient pml_c2_Re(dim, detJ_JT_J_inv_Re, pml);
    PMLDiagMatrixCoefficient pml_c2_Im(dim, detJ_JT_J_inv_Im, pml);
    mfem::ScalarVectorProductCoefficient c2_Re(omeg, pml_c2_Re);
    mfem::ScalarVectorProductCoefficient c2_Im(omeg, pml_c2_Im);
    mfem::VectorRestrictedCoefficient restr_c2_Re(c2_Re, attrPML);
    mfem::VectorRestrictedCoefficient restr_c2_Im(c2_Im, attrPML);

    // Integrators inside the PML region
    a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(restr_c1_Re),
        new mfem::CurlCurlIntegrator(restr_c1_Im));
    a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(restr_c2_Re),
        new mfem::VectorFEMassIntegrator(restr_c2_Im));

    // Assemble the bilinear form and the corresponding linear system
    if (pa) { a.SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL); }
    a.Assemble(0);

    mfem::OperatorPtr A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // 14. Solving and preconditioning

    if (pa || !umf_solver)
    {
        mfem::ConstantCoefficient absomeg(pow2(omega) * epsilon);
        mfem::RestrictedCoefficient restr_absomeg(absomeg, attr);

        mfem::BilinearForm prec(fespace);
        prec.AddDomainIntegrator(new mfem::CurlCurlIntegrator(restr_muinv));
        prec.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(restr_absomeg));

        PMLDiagMatrixCoefficient pml_c1_abs(cdim, detJ_inv_JT_J_abs, pml);
        mfem::ScalarVectorProductCoefficient c1_abs(muinv, pml_c1_abs);
        mfem::VectorRestrictedCoefficient restr_c1_abs(c1_abs, attrPML);

        PMLDiagMatrixCoefficient pml_c2_abs(dim, detJ_JT_J_inv_abs, pml);
        mfem::ScalarVectorProductCoefficient c2_abs(absomeg, pml_c2_abs);
        mfem::VectorRestrictedCoefficient restr_c2_abs(c2_abs, attrPML);

        prec.AddDomainIntegrator(new mfem::CurlCurlIntegrator(restr_c1_abs));
        prec.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(restr_c2_abs));

        if (pa) { prec.SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL); }
        prec.Assemble();

        // Define and apply a GMRES solver for AU=B with a block diagonal
        //      preconditioner based on the Gauss-Seidel or Jacobi sparse smoother.
        mfem::Array<int> offsets(3);
        offsets[0] = 0;
        offsets[1] = fespace->GetTrueVSize();
        offsets[2] = fespace->GetTrueVSize();
        offsets.PartialSum();

        std::unique_ptr<mfem::Operator> pc_r;
        std::unique_ptr<mfem::Operator> pc_i;
        real_t s = (conv == mfem::ComplexOperator::HERMITIAN) ? -1_r : 1_r;
        if (pa)
        {
            // Jacobi Smoother
            pc_r.reset(new mfem::OperatorJacobiSmoother(prec, ess_tdof_list));
            pc_i.reset(new mfem::ScaledOperator(pc_r.get(), s));
        }
        else
        {
            mfem::OperatorPtr PCOpAh;
            prec.SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
            prec.FormSystemMatrix(ess_tdof_list, PCOpAh);

            // Gauss-Seidel Smoother
            pc_r.reset(new mfem::GSSmoother(*PCOpAh.As<mfem::SparseMatrix>()));
            pc_i.reset(new mfem::ScaledOperator(pc_r.get(), s));
        }

        mfem::BlockDiagonalPreconditioner BlockDP(offsets);
        BlockDP.SetDiagonalBlock(0, pc_r.get());
        BlockDP.SetDiagonalBlock(1, pc_i.get());

        mfem::GMRESSolver gmres;
        gmres.SetPrintLevel(1);
        gmres.SetKDim(200);
        gmres.SetMaxIter(pa ? 5000 : 2000);
        gmres.SetRelTol(1e-5);
        gmres.SetAbsTol(0.0);
        gmres.SetOperator(*A);
        gmres.SetPreconditioner(BlockDP);
        gmres.Mult(B, X);
    }

    // Recover the solution as a finite element grid function.
    a.RecoverFEMSolution(X, b, x);

    // 16. Save the refined mesh and the solution, viewable using GLViS.
    std::ofstream sol_r_ofs("sol_r.gf");
    std::ofstream sol_i_ofs("sol_i.gf");
    sol_r_ofs.precision(8);
    sol_i_ofs.precision(8);
    x.real().Save(sol_r_ofs);
    x.imag().Save(sol_i_ofs);

    if (visualization)
    {
        // Define visualization keys for GLVis (see GLVis documentation)
        std::string keys;
        keys = (dim == 3) ? "keys macF\n" : keys = "keys amrRljcUUuu\n";

        char vishost[] = "localhost";
        int visport = 19916;

        mfem::socketstream sol_sock_re(vishost, visport);
        sol_sock_re.precision(8);
        sol_sock_re << "solution\n"
            << *mesh << x.real() << keys
            << "window_title 'Solution real part'" << std::flush;

        mfem::socketstream sol_sock_im(vishost, visport);
        sol_sock_im.precision(8);
        sol_sock_im << "solution\n"
            << *mesh << x.imag() << keys
            << "window_title 'Solution imag part'" << std::flush;

        mfem::GridFunction x_t(fespace);
        x_t = x.real();
        mfem::socketstream sol_sock(vishost, visport);
        sol_sock.precision(8);
        sol_sock << "solution\n"
            << *mesh << x_t << keys << "autoscale off\n"
            << "window_title 'Harmonic Solution (t = 0.0 T)'"
            << "pause\n" << std::flush;
        std::cout << "GLVis visualization paused."
            << " Press space (in the GLVis window) to resume it.\n";
        int num_frames = 32;
        int i = 0;
        while (sol_sock)
        {
            real_t t = (real_t)(i % num_frames) / num_frames;
            std::ostringstream oss;
            oss << "Harmonic Solution (t = " << t << " T)";

            add(cos(real_t(2.0 * M_PI) * t), x.real(),
                sin(real_t(2.0 * M_PI) * t), x.imag(), x_t);
            sol_sock << "solution\n"
                << *mesh << x_t
                << "window_title '" << oss.str() << "'" << std::flush;
            i++;
        }
    }
}
void ForwardSolve::set_inverse_system()
{
    /*
    * INVERSE MATRICES SETUP
    */

    // ****************************************************************************************
    // I. e vectors: 1. Create elementwise e field vectors
    std::vector<mfem::Vector> e_real_locals;
    std::vector<mfem::Vector> e_imag_locals;

    for (int i = 0; i < fespace_x->GetNE(); i++)
    {
        mfem::Array<int> dofs;
        fespace_x->GetElementDofs(i, dofs); // Get DOFs for this element

        // Create a local element vector
        mfem::Vector e_r(dofs.Size());
        mfem::Vector e_i(dofs.Size());

        x.real().GetSubVector(dofs, e_r);
        x.imag().GetSubVector(dofs, e_i);

        // Store the local element vector in the container
        e_real_locals.push_back(e_r);
        e_imag_locals.push_back(e_i);
    }

    // ****************************************************************************************
    // II. M matrix: 1. Obtaining the Mass matrix related to the time derivative term
    mfem::SesquilinearForm M(fespace, conv);
    M.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(restr_sq_omega_eps), NULL);
    M.Assemble(0);
    mfem::SparseMatrix& M_real = M.real().SpMat(); // Real part of the mass matrix
    mfem::SparseMatrix& M_imag = M.imag().SpMat(); // Imaginary part of the mass matrix
    M_real.Finalize();
    M_imag.Finalize();
    //std::cout << "Size of M_real: " << M_real.Height() << " x " << M_real.Width() << std::endl;
 //   M_real.Print(std::cout);
 //   std::cout << "Size of M_imag: " << M_imag.Height() << " x " << M_imag.Width() << std::endl;
 //   M_imag.Print(std::cout);

    // ****************************************************************************************
    // II. M matrix: 2. Using Ms to obtain the elementwise M matrices

    // Define a container to store local elementwise matrices
    std::vector<mfem::DenseMatrix> M_real_locals;
    std::vector<mfem::DenseMatrix> M_imag_locals;

    // Loop over all elements in the mesh
    for (int i = 0; i < mesh->GetNE(); i++)
    {
        // Get the finite element associated with the element
        const mfem::FiniteElement* fe = fespace->GetFE(i);

        // Get the degrees of freedom for this element
        mfem::Array<int> dofs;
        fespace->GetElementDofs(i, dofs);

        // Create a local element mass matrix
        mfem::DenseMatrix m_r(dofs.Size(), dofs.Size());
        mfem::DenseMatrix m_i(dofs.Size(), dofs.Size());

        // Might need to add/subtract boundary values between elements
        M_real.GetSubMatrix(dofs, dofs, m_r);
        M_imag.GetSubMatrix(dofs, dofs, m_i);

        // Store the local element mass matrices
        M_real_locals.push_back(m_r);
        M_imag_locals.push_back(m_i);
    }

    // ****************************************************************************************
    // III. Setting Up Imaging Region DOFs
    double sim_region_length = 2.0;

    // Store DOFs of all elements within the imaging region
    std::vector<int> region_elements;
    std::vector<mfem::Array<int>> region_dofs;

    for (int i = 0; i < mesh->GetNE(); i++)
    {
        mfem::ElementTransformation* trans = mesh->GetElementTransformation(i);
        mfem::Vector centroid(mesh->Dimension());
        trans->Transform(mfem::Geometries.GetCenter(mesh->GetElementBaseGeometry(i)), centroid);

        if (IsWithinSquare(centroid, circle_center, sim_region_length / 2))
        {
            region_elements.push_back(i);
            mfem::Array<int> dofs;
            fespace->GetElementDofs(i, dofs);
            region_dofs.push_back(dofs);
        }
    }

    // ****************************************************************************************
    // IV. Setting up Locations of Elements of Receivers and Sources

    std::vector<int> receivers_elements;
    std::vector<mfem::Array<int>> receivers_dofs;
    int source_element{ -1 };
    mfem::Array<int> source_dofs;

    for (int i = 0; i < system.num_antennas; i++)
    {
        if (i == system.iteration)
        {
            mfem::Vector a(dim);
            a(0) = system.location_source.x.at(i);
            a(1) = system.location_source.y.at(i);
            double min_dist = std::numeric_limits<double>::max();

            for (int j = 0; j < mesh->GetNE(); j++)
            {
                mfem::ElementTransformation* trans = mesh->GetElementTransformation(j);
                mfem::Vector center(dim);
                trans->Transform(mfem::Geometries.GetCenter(mesh->GetElementBaseGeometry(j)), center);

                double dist = center.DistanceTo(a);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    source_element = j;
                }
            }
            fespace->GetElementDofs(source_element, source_dofs);
        }
        else
        {
            mfem::Vector a(dim);
            a(0) = system.location_receivers.x.at(i - 1);
            a(1) = system.location_receivers.y.at(i - 1);
            double min_dist = std::numeric_limits<double>::max();
            int closest_elem = -1;

            for (int j = 0; j < mesh->GetNE(); j++)
            {
                mfem::ElementTransformation* trans = mesh->GetElementTransformation(j);
                mfem::Vector center(dim);
                trans->Transform(mfem::Geometries.GetCenter(mesh->GetElementBaseGeometry(j)), center);

                double dist = center.DistanceTo(a);
                if (dist < min_dist)
                {
                    min_dist = dist;
                    closest_elem = j;
                }
            }

            receivers_elements.push_back(closest_elem);

            // Get DOFs for the closest element
            mfem::Array<int> dofs;
            fespace->GetElementDofs(closest_elem, dofs);
            receivers_dofs.push_back(dofs);
        }
    }

    // ****************************************************************************************
    // V. e_j vector: 1. Extracting receiver and source location specific e field values

    std::vector<mfem::Vector> e_real_j;
    std::vector<mfem::Vector> e_imag_j;
    std::vector<mfem::Vector> e_real_tx;
    std::vector<mfem::Vector> e_imag_tx;

    std::vector<mfem::Vector> e_j_ri;

    for (int i = 0; i < system.num_antennas; i++)
    {
        if (i == system.iteration)
        {
            mfem::Vector e_r(source_dofs.Size());
            mfem::Vector e_i(source_dofs.Size());

            e_real_locals.at(source_element).GetSubVector(source_dofs, e_r);
            e_imag_locals.at(source_element).GetSubVector(source_dofs, e_i);

            e_real_tx.push_back(e_r);
            e_imag_tx.push_back(e_i);

        }
        else
        {
            mfem::Vector e_r(receivers_dofs.at(i - 1).Size());
            mfem::Vector e_i(receivers_dofs.at(i - 1).Size());

            e_real_locals.at(receivers_elements.at(i - 1)).GetSubVector(receivers_dofs.at(i - 1), e_r);
            e_imag_locals.at(receivers_elements.at(i - 1)).GetSubVector(receivers_dofs.at(i - 1), e_i);

            e_real_j.push_back(e_r);
            e_imag_j.push_back(e_i);
            e_j_ri.push_back(e_r);
            e_j_ri.push_back(e_i);
        }
    }

    // ****************************************************************************************
    // VI. e_i vector: 1. Extracting imaging region specific e field values

    std::vector<mfem::Vector> e_real_i;
    std::vector<mfem::Vector> e_imag_i;
    std::vector<mfem::Vector> e_i_ri;

    for (int i = 0; i < region_elements.size(); i++)
    {
        mfem::Vector e_r(region_dofs.at(i).Size());
        mfem::Vector e_i(region_dofs.at(i).Size());

        e_real_locals.at(region_elements.at(i)).GetSubVector(region_dofs.at(i), e_r);
        e_imag_locals.at(region_elements.at(i)).GetSubVector(region_dofs.at(i), e_i);

        e_real_i.push_back(e_r);
        e_imag_i.push_back(e_i);
        e_i_ri.push_back(e_r);
        e_i_ri.push_back(e_i);
    }



    // ****************************************************************************************
    // VII. M matrix: 3. Extracting imaging region specific M matrices

    std::vector<mfem::DenseMatrix> M_real_i;
    std::vector<mfem::DenseMatrix> M_imag_i;

    for (int i = 0; i < region_elements.size(); i++)
    {
        mfem::DenseMatrix m_r(region_dofs.at(i).Size(), region_dofs.at(i).Size());
        mfem::DenseMatrix m_i(region_dofs.at(i).Size(), region_dofs.at(i).Size());

        M_real.GetSubMatrix(region_dofs.at(i), region_dofs.at(i), m_r);
        M_imag.GetSubMatrix(region_dofs.at(i), region_dofs.at(i), m_i);

        M_real_i.push_back(m_r);
        M_imag_i.push_back(m_i);
    }

    // ****************************************************************************************
    // III. E_ref vector: 1. Reading the reference data
    mfem::GridFunction E_ref_real(fespace_x), E_ref_imag(fespace_x);
    std::ifstream in("ex25-sol_r.gf");
    E_ref_real.Load(in);
    in.close();
    in.open("ex25-sol_i.gf");
    E_ref_imag.Load(in);
    in.close();

    // ****************************************************************************************
    // IV. E_ref vector: 2. Extracting transmitter and receiver location specific e field values from Part I

    // Define a container to store the receiver location specific e field values

    for (int i = 0; i < system.num_antennas - 1; i++)
    {
        mfem::Vector e_r(receivers_dofs.at(i).Size());
        mfem::Vector e_i(receivers_dofs.at(i).Size());

        for (int j = 0; j < receivers_dofs.at(i).Size(); j++)
        {
            e_r(j) = E_ref_real(receivers_dofs.at(i)[j]);
            e_i(j) = E_ref_imag(receivers_dofs.at(i)[j]);
        }
        e_ref_real_rx.push_back(e_r);
        e_ref_imag_rx.push_back(e_i);
    }
    // V. Setting the A matrix for inversion: A_{m,n} = \sum_{j=1}^{N} e_{m,j}^T M_j e_{n,j}

    // Initialising required matrices and vectors
    size_t num_dofs = region_dofs.at(0).Size();

    mfem::Vector e_j_n(2 * num_dofs);
    mfem::DenseMatrix e_i_n(2 * num_dofs, 2 * num_dofs);
    mfem::DenseMatrix M_n(2 * num_dofs, 2 * num_dofs);

    mfem::Vector A_real_n(num_dofs);
    mfem::Vector A_imag_n(num_dofs);

    size_t num_rows = receivers_elements.size();
    size_t num_cols = region_elements.size() * num_dofs;

    mfem::SparseMatrix A_real(num_rows, num_cols);
    mfem::SparseMatrix A_imag(num_rows, num_cols);

    A_inv = mfem::SparseMatrix(2 * num_rows, 2 * num_cols);

    // V. Setting the elementwise A matrix
    for (int m = 0; m < receivers_elements.size(); m++)
    {
        // e_j vector for the m-th receiver
        for (int i = 0; i < 2 * num_dofs; i++)
        {
            e_j_n(i) = e_j_ri.at(m)(i);
        }

        // e_i and M matrices for the region
        for (int n = 0; n < region_elements.size(); n++)
        {
            for (int j = 0; j < num_dofs; j++)
            {
                e_i_n(0, j) = e_real_i.at(n)[j];                // Top-left
                e_i_n(0, j + num_dofs) = e_imag_i.at(n)[j];     // Top-right
                e_i_n(1, j) = -e_imag_i.at(n)[j];               // Bottom-left
                e_i_n(1, j + num_dofs) = e_real_i.at(n)[j];     // Bottom-right

                M_n(0, j) = M_real_i.at(n)(j, 0);                // Top-left
                M_n(0, j + num_dofs) = M_imag_i.at(n)(j, 0);     // Top-right
                M_n(1, j) = -M_imag_i.at(n)(j, 0);               // Bottom-left
                M_n(1, j + num_dofs) = M_real_i.at(n)(j, 0);     // Bottom-right
            }

            // A_{m,n} = e_{m,j}^T M_j e_{n,j}
            mfem::Vector temp(2 * num_dofs);
            M_n.Mult(e_j_n, temp);

            mfem::Vector A_n(2 * num_dofs);
            e_i_n.Mult(temp, A_n);

            // Separate the elementwise A matrix
            for (int i = 0; i < num_dofs; i++)
            {
                A_real_n(i) = A_n(i);
                A_imag_n(i) = A_n(i + num_dofs);
            }

            // Store the elementwise A matrix
            A_real.SetRow(m, region_dofs.at(n), A_real_n);
            A_imag.SetRow(m, region_dofs.at(n), A_imag_n);
        }
    }

    A_real.Finalize();
    A_imag.Finalize();

    mfem::SparseMatrix A_real_minus(num_rows, num_cols);
    A_real_minus = 0.0;
    A_real_minus.Add(-1, A_real);

    // Assembling the final A matrix by: A = [A_real, A_imag; A_imag, -A_real]
    for (int i = 0; i < num_rows; i++)
    {
        mfem::Vector A_real_row(num_cols);
        mfem::Vector A_imag_row(num_cols);
        mfem::Vector A_real_minus_row(num_cols);

        A_real.GetRow(i, region_dofs.at(0), A_real_row);
        A_imag.GetRow(i, region_dofs.at(0), A_imag_row);
        A_real_minus.GetRow(i, region_dofs.at(0), A_real_minus_row);

        A_inv.SetRow(i, region_dofs.at(0), A_real_row);
        A_inv.SetRow(i + num_rows, region_dofs.at(0), A_imag_row);

        A_inv.SetRow(i, region_dofs.at(region_dofs.size() - 1), A_imag_row);
        A_inv.SetRow(i + num_rows, region_dofs.at(region_dofs.size() - 1), A_real_minus_row);
    }

    A_inv.Finalize();

    // VI. Setting the b vector for inversion
    b_inv = mfem::Vector(2 * num_rows);

    for (int i = 0; i < num_rows; i++)
    {
        b_inv(i) = e_ref_real_rx.at(i) - e_real_j.at(i);
        b_inv(i + num_rows) = e_ref_imag_rx.at(i) - e_imag_j.at(i);
    }
}
double ForwardSolve::evaluate_residual()
{
    // Residual: r = b/sqrt(E_ref_real^2+E_ref_imag^2)
	mfem::Vector r = b_inv;

	for (int i = 0; i < r.Size(); i++)
	{
		for (int j = 0; j < e_ref_real_rx.size(); j++)
		{
            r(i) = r(i) / sqrt(pow2(e_ref_real_rx.at(i)(j)) + pow2(e_ref_imag_rx.at(i)(j)));
		}
	}

	return r.Norml2();
}

mfem::Vector ForwardSolve::evaluate_residual_derivative()
{
	// Derivative of the residual: dr = (A^T * b) / sqrt(E_ref_real^2+E_ref_imag^2)

    // Incomplete
	mfem::Vector dr(2 * A_inv.Height());

	for (int i = 0; i < dr.Size(); i++)
	{
		for (int j = 0; j < e_ref_real_rx.size(); j++)
		{
			dr(i) = dr(i) / sqrt(pow2(e_ref_real_rx.at(i)(j)) + pow2(e_ref_imag_rx.at(i)(j)));
		}
	}
}


// MFEM Forward Solve: Constituent Functions
bool ForwardSolve::IsWithinSquare(const mfem::Vector& point, const mfem::Vector& center, double half_size)
{
    return (point(0) >= center(0) - half_size) &&
        (point(0) <= center(0) + half_size) &&
        (point(1) >= center(1) - half_size) &&
        (point(1) <= center(1) + half_size);
}

void ForwardSolve::eval_source_receivers_locations(
    const size_t num_antennas,
    const double radius,
    const size_t iteration,
    Points2D<double>& location_source,
    Points2D<double>& location_receivers)
{
    for (size_t i = 0; i < num_antennas; i++)
    {
        double angle = 2 * M_PI * i / num_antennas;
        if (i == iteration)
        {
            system.location_source.x.push_back(radius * cos(angle));
            system.location_source.y.push_back(radius * sin(angle));
        }
        else
        {
            system.location_receivers.x.push_back(radius * cos(angle));
            system.location_receivers.y.push_back(radius * sin(angle));
        }
    }
}

void ForwardSolve::source(const mfem::Vector& x, mfem::Vector& f)
{
    // INVERSE SOURCE RECEIVER LOCATIONS AND FUNCTION SETUP
    eval_source_receivers_locations(system.num_antennas, system.radius, system.iteration, system.location_source, system.location_receivers);

    mfem::Vector center(dim);
    real_t r = 0.0;
    center(0) = system.location_source.x.at(0);
    center(1) = system.location_source.y.at(0);
    for (int i = 0; i < dim; ++i)
    {
        center(i) = 0.5_r * (comp_domain_bdr(i, 0) + comp_domain_bdr(i, 1));
        r += pow2(x[i] - center[i]);
    }
    // Implement the delta pulse as a Gaussian source
    real_t n = 5_r * omega * sqrt(epsilon * mu) / real_t(M_PI);
    real_t coeff = pow2(n) / real_t(M_PI);
    real_t alpha = -pow2(n) * r;
    f = 0.0;
    f[0] = coeff * exp(alpha);
}

void maxwell_solution(const mfem::Vector& x, std::vector<std::complex<real_t>>& E)
{
    // Initialize
    for (int i = 0; i < dim; ++i)
    {
        E[i] = 0.0;
    }

    constexpr std::complex<real_t> zi = std::complex<real_t>(0., 1.);
    real_t k = omega * sqrt(epsilon * mu);
}

void E_data_Re(const mfem::Vector& x, mfem::Vector& E)
{
    std::vector<std::complex<real_t>> Eval(E.Size());
    for (int i = 0; i < dim; ++i)
    {
        E[i] = Eval[i].real();
    }
}

void E_data_Im(const mfem::Vector& x, mfem::Vector& E)
{
    std::vector<std::complex<real_t>> Eval(E.Size());
    for (int i = 0; i < dim; ++i)
    {
        E[i] = Eval[i].imag();
    }
}

void E_bdr_data_Re(const mfem::Vector& x, mfem::Vector& E)
{
    E = 0.0;
    bool in_pml = false;

    for (int i = 0; i < dim; ++i)
    {
        // check if in PML
        if (x(i) - comp_domain_bdr(i, 0) < 0_r ||
            x(i) - comp_domain_bdr(i, 1) > 0_r)
        {
            in_pml = true;
            break;
        }
    }
    if (!in_pml)
    {
        std::vector<std::complex<real_t>> Eval(E.Size());
        maxwell_solution(x, Eval);
        for (int i = 0; i < dim; ++i)
        {
            E[i] = Eval[i].real();
        }
    }
}

// Define bdr_data solution
void E_bdr_data_Im(const mfem::Vector& x, mfem::Vector& E)
{
    E = 0.0;
    bool in_pml = false;

    for (int i = 0; i < dim; ++i)
    {
        // check if in PML
        if (x(i) - comp_domain_bdr(i, 0) < 0_r ||
            x(i) - comp_domain_bdr(i, 1) > 0_r)
        {
            in_pml = true;
            break;
        }
    }
    if (!in_pml)
    {
        std::vector<std::complex<real_t>> Eval(E.Size());
        maxwell_solution(x, Eval);
        for (int i = 0; i < dim; ++i)
        {
            E[i] = Eval[i].imag();
        }
    }
}

void detJ_JT_J_inv_Re(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det(1.0, 0.0);
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    for (int i = 0; i < dim; ++i)
    {
        D(i) = (det / pow2(dxs[i])).real();
    }
}

void detJ_JT_J_inv_Im(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det = 1.0;
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    for (int i = 0; i < dim; ++i)
    {
        D(i) = (det / pow2(dxs[i])).imag();
    }
}

void detJ_JT_J_inv_abs(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det = 1.0;
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    for (int i = 0; i < dim; ++i)
    {
        D(i) = abs(det / pow2(dxs[i]));
    }
}

void detJ_inv_JT_J_Re(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det(1.0, 0.0);
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    // in the 2D case the coefficient is scalar 1/det(J)
    if (dim == 2)
    {
        D = (1_r / det).real();
    }
    else
    {
        for (int i = 0; i < dim; ++i)
        {
            D(i) = (pow2(dxs[i]) / det).real();
        }
    }
}

void detJ_inv_JT_J_Im(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det = 1.0;
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    if (dim == 2)
    {
        D = (1_r / det).imag();
    }
    else
    {
        for (int i = 0; i < dim; ++i)
        {
            D(i) = (pow2(dxs[i]) / det).imag();
        }
    }
}

void detJ_inv_JT_J_abs(const mfem::Vector& x, PML* pml, mfem::Vector& D)
{
    std::vector<std::complex<real_t>> dxs(dim);
    std::complex<real_t> det = 1.0;
    pml->StretchFunction(x, dxs);

    for (int i = 0; i < dim; ++i)
    {
        det *= dxs[i];
    }

    if (dim == 2)
    {
        D = abs(1_r / det);
    }
    else
    {
        for (int i = 0; i < dim; ++i)
        {
            D(i) = abs(pow2(dxs[i]) / det);
        }
    }
}

PML::PML(mfem::Mesh* mesh_, mfem::Array2D<real_t> length_)
    : mesh(mesh_), length(length_)
{
    dim = mesh->Dimension();
    SetBoundaries();
}

void PML::SetBoundaries()
{
    comp_dom_bdr.SetSize(dim, 2);
    dom_bdr.SetSize(dim, 2);
    mfem::Vector pmin, pmax;
    mesh->GetBoundingBox(pmin, pmax);
    for (int i = 0; i < dim; i++)
    {
        dom_bdr(i, 0) = pmin(i);
        dom_bdr(i, 1) = pmax(i);
        comp_dom_bdr(i, 0) = dom_bdr(i, 0) + length(i, 0);
        comp_dom_bdr(i, 1) = dom_bdr(i, 1) - length(i, 1);
    }
}

void PML::SetAttributes(mfem::Mesh* mesh_)
{
    // Initialize bdr attributes
    for (int i = 0; i < mesh_->GetNBE(); ++i)
    {
        mesh_->GetBdrElement(i)->SetAttribute(i + 1);
    }

    int nrelem = mesh_->GetNE();

    elems.SetSize(nrelem);

    // Loop through the elements and identify which of them are in the PML
    for (int i = 0; i < nrelem; ++i)
    {
        elems[i] = 1;
        bool in_pml = false;
        mfem::Element* el = mesh_->GetElement(i);
        mfem::Array<int> vertices;

        // Initialize attribute
        el->SetAttribute(1);
        el->GetVertices(vertices);
        int nrvert = vertices.Size();

        // Check if any vertex is in the PML
        for (int iv = 0; iv < nrvert; ++iv)
        {
            int vert_idx = vertices[iv];
            real_t* coords = mesh_->GetVertex(vert_idx);
            for (int comp = 0; comp < dim; ++comp)
            {
                if (coords[comp] > comp_dom_bdr(comp, 1) ||
                    coords[comp] < comp_dom_bdr(comp, 0))
                {
                    in_pml = true;
                    break;
                }
            }
        }
        if (in_pml)
        {
            elems[i] = 0;
            el->SetAttribute(2);
        }
    }
    mesh_->SetAttributes();
}

void PML::StretchFunction(const mfem::Vector& x,
    std::vector<std::complex<real_t>>& dxs)
{
    constexpr std::complex<real_t> zi = std::complex<real_t>(0., 1.);

    real_t n = 2.0;
    real_t c = 5.0;
    real_t coeff;
    real_t k = omega * sqrt(epsilon * mu);

    // Stretch in each direction independently
    for (int i = 0; i < dim; ++i)
    {
        dxs[i] = 1.0;
        if (x(i) >= comp_domain_bdr(i, 1))
        {
            coeff = n * c / k / pow(length(i, 1), n);
            dxs[i] = 1_r + zi * coeff *
                abs(pow(x(i) - comp_domain_bdr(i, 1), n - 1_r));
        }
        if (x(i) <= comp_domain_bdr(i, 0))
        {
            coeff = n * c / k / pow(length(i, 0), n);
            dxs[i] = 1_r + zi * coeff *
                abs(pow(x(i) - comp_domain_bdr(i, 0), n - 1_r));
        }
    }
}