#pragma once
#include "DefaultMesh.h"

// Class for adaptive mesh generation inherited from DefaultMesh
class AdaptiveMesh : public DefaultMesh
{
public:

	float scaling_factor;
	std::string refine_type;

public:
	AdaptiveMesh(MWInvSystemSetup& mw_sys, std::string mesh_name, bool display_mesh, bool print_data);
	~AdaptiveMesh();
	AdaptiveMesh& operator=(const AdaptiveMesh& mesh);
	// Functions overloaded from DefaultMesh
	// Function to generate mesh, given user inputs of the system and options to display mesh and print data
	void generate_mesh(std::string mesh_name, bool display_mesh, bool print_data);
	// Function to print information from the mesh
	void print_mesh_information();
};

