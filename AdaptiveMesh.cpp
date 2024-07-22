#include "AdaptiveMesh.h"

AdaptiveMesh::AdaptiveMesh(MWInvSystemSetup& mw_sys, std::string mesh_name, bool display_mesh, bool print_data)
	:
scaling_factor(1.0), refine_type("uniform")
{}
AdaptiveMesh::~AdaptiveMesh()
{}
void AdaptiveMesh::generate_mesh()
{}
void AdaptiveMesh::print_mesh_information()
{}
AdaptiveMesh& AdaptiveMesh::operator=(const AdaptiveMesh& mesh)
{}