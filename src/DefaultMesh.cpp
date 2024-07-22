#include "DefaultMesh.h"

DefaultMesh::DefaultMesh(
	MWInvSystemSetup & mw_sys, 
	std::string mesh_name, 
	bool display_mesh, 
	bool print_data) : 
	mw_sys(mw_sys),
	mesh_name(mesh_name),
	display_mesh(display_mesh),
	print_data(print_data)
{}

DefaultMesh::~DefaultMesh() {}

DefaultMesh& DefaultMesh::operator=(const DefaultMesh& mesh)
{
	if (this == &mesh)
	{
		return *this;
	}

	mw_sys = mesh.mw_sys;
	mesh_name = mesh.mesh_name;
	display_mesh = mesh.display_mesh;
	print_data = mesh.print_data;

	return *this;
}

void DefaultMesh::print_mesh_information()
{
	// Extract node data
	std::vector<std::size_t> nodeTags;
	std::vector<double> nodeCoords;
	std::vector<double> parametricCoords;

	gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

	// Print node data
	std::cout << "Number of Nodes:" << nodeTags.size() << std::endl;

	// Extract element data
	std::vector<int> elementTypes;
	std::vector<std::vector<std::size_t>> elementTags, elementNodeTags;

	gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
	std::cout << "Number of Elements:" << elementTags.size() << std::endl;
}

void DefaultMesh::generate_mesh(std::string mesh_name, bool display_mesh, bool print_data)
{
    gmsh::initialize();
    gmsh::model::add(mesh_name);

	double lc = 1e-2; // Evaluate this using Wave propagation and adaptive meshing

    /*Adding All Relevant Points*/

	// Adding source and receiver points
	for (size_t i = 0; i < mw_sys.num_antennas; i++)
    {
		if (i == mw_sys.iteration) { gmsh::model::geo::addPoint(mw_sys.location_source.x.at(i), mw_sys.location_source.y.at(i), 0, lc, i + 1); }
		else { gmsh::model::geo::addPoint(mw_sys.location_receivers.x.at(i), mw_sys.location_receivers.y.at(i), 0, lc, i + 1); }
	}

	// Adding the maximum simulation region boundary points
	gmsh::model::geo::addPoint(-mw_sys.max_region_length / 2, -mw_sys.max_region_width / 2, 0, lc, mw_sys.num_antennas + 1);
	gmsh::model::geo::addPoint(mw_sys.max_region_length / 2, -mw_sys.max_region_width / 2, 0, lc, mw_sys.num_antennas + 2);
	gmsh::model::geo::addPoint(mw_sys.max_region_length / 2, mw_sys.max_region_width / 2, 0, lc, mw_sys.num_antennas + 3);
	gmsh::model::geo::addPoint(-mw_sys.max_region_length / 2, mw_sys.max_region_width / 2, 0, lc, mw_sys.num_antennas + 4);

	/*Adding All Relevant Lines and Curves*/

	// Adding maximum region boundary lines
	gmsh::model::geo::addLine(mw_sys.num_antennas + 1, mw_sys.num_antennas + 2, 1);
	gmsh::model::geo::addLine(mw_sys.num_antennas + 2, mw_sys.num_antennas + 3, 2);
	gmsh::model::geo::addLine(mw_sys.num_antennas + 3, mw_sys.num_antennas + 4, 3);
	gmsh::model::geo::addLine(mw_sys.num_antennas + 4, mw_sys.num_antennas + 1, 4);

    /*Adding Relevant Surfaces*/

	// Adding the maximum region boundary surface
	gmsh::model::geo::addCurveLoop({ 1, 2, 3, 4 }, 1);
	gmsh::model::geo::addPlaneSurface({ 1 }, 1);

    gmsh::model::geo::synchronize();

    /*Defining the Physical Groups*/

	// Defining All Source Points as One Physical Group
	std::vector<int> antenna_points;
	for (size_t i = 1; i <= mw_sys.num_antennas; i++)
	{
		antenna_points.push_back(i);
	}

	gmsh::model::addPhysicalGroup(0, antenna_points, 1);

	// Defining the Maximum Region Boundary as One Physical Group
	gmsh::model::addPhysicalGroup(2, { 1 }, 3);

	// Generate the mesh
    gmsh::model::mesh::generate(2);

	// Save the mesh
	std::string mesh_file_name = mesh_name + ".msh";
	gmsh::write(mesh_file_name);

	if (display_mesh)
	{ gmsh::fltk::run(); }

	if (print_data)
	{ print_mesh_information();}

    gmsh::finalize();
}