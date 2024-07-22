#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "MFEM/mfem.hpp"

#include "Points2D.h"

class MWInvSystemSetup
{
public:
	// Simulation region maximmum extents
	const double max_region_length;
	const double max_region_width;
	// Region inside the max region where the reconstruction will be done
	const double sim_region_length;
	const double sim_region_width;

	//Source and receiver
	// Number of sources and receivers
	const size_t  num_antennas;
	const size_t iteration;

	// Distance - Radius from the centre of region 
	// to where sources and receivers are placed
	const double radius;

	// Locations of sources and receivers
	Points2D<double> location_source;
	Points2D<double> location_receivers;
	
	// Source waveform parameters
	std::vector<double> source_intensities;
	std::vector<double> source_frequencies;
	std::string source_type;

	// Region permittivity parameters
	static const double epsilon_0;
	std::string object_type;	// Based on tests I will do
	double object_epsilon_r;
	double object_sigma;
	std::vector<double> object_bounds;	// Will be radius for circle, length and width for rectangle, etc.

public:
	MWInvSystemSetup(); //Stick to default constructor
	~MWInvSystemSetup();

	// Function to take user inputs and return the system setup
	void input_system_setup();
	
	// Function to evaluate the source and receivers locations
	void evaluate_source_receivers_locations(const size_t num_antennas, const double radius, const size_t iteration, Points2D<double>& location_source, Points2D<double>& location_receivers);

	// Function to print the system information
	void print_system_information();
	
	// Overloaded constructor to set the system setup based on the inputs
	MWInvSystemSetup(
		double max_region_length,
		double max_region_width,
		double sim_region_length,
		double sim_region_width,
		size_t num_antennas,
		double radius,
		size_t iteration,
		std::string source_type,
		std::vector<double> source_intensities,
		std::vector<double> source_frequencies,
		std::string object_type,
		double object_epsilon_r,
		double object_sigma,
		std::vector<double> object_bounds);

	// Copy constructor
	MWInvSystemSetup(const MWInvSystemSetup& mw_sys);

	// Assignment operator
	MWInvSystemSetup& operator=(const MWInvSystemSetup& mw_sys);

};