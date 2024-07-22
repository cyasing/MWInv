#pragma once

#include "MWInvSystemSetup.h"
#include "gmsh/gmsh.h"

// Generates Mesh for the system based on user inputs
class DefaultMesh
{
private:
	MWInvSystemSetup& mw_sys;
	std::string mesh_name;
	bool display_mesh;
	bool print_data;
public:
	DefaultMesh(MWInvSystemSetup& mw_sys, std::string mesh_name, bool display_mesh, bool print_data);
	~DefaultMesh();
	// operator=
	DefaultMesh& operator=(const DefaultMesh& mesh);
	
	// Function to generate mesh, given user inputs of the system and options to display mesh and print data
	void generate_mesh(std::string mesh_name, bool display_mesh, bool print_data);
	// Function to print information from the mesh
	void print_mesh_information();
};