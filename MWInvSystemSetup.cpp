#include "MWInvSystemSetup.h"

#define M_PI 3.14159265358979323846
const double MWInvSystemSetup::epsilon_0 = 8.854187817e-12;

MWInvSystemSetup::MWInvSystemSetup() :
	max_region_length(4.0), 
	max_region_width(4.0), 
	sim_region_length(2.0), 
	sim_region_width(2.0), 
	num_antennas(8), 
	radius(1.8),
	iteration(0),
	source_frequencies({ 1.0, 2.0, 3.0, 4.0 }),
	source_type("sine"),
	source_intensities({ 1.0, 1.0, 1.0, 1.0 }),
	object_type("circle"), 
	object_epsilon_r(4.0), 
	object_sigma(0.0), 
	object_bounds({ 0.5 })
{
	evaluate_source_receivers_locations(num_antennas,radius,iteration, location_source, location_receivers);
}

MWInvSystemSetup::MWInvSystemSetup(
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
	std::vector<double> object_bounds) :

	max_region_length(max_region_length),
	max_region_width(max_region_width),
	sim_region_length(sim_region_length),
	sim_region_width(sim_region_width),
	num_antennas(num_antennas),
	radius(radius),
	iteration(iteration),
	object_type(object_type),
	object_epsilon_r(object_epsilon_r),
	object_sigma(object_sigma),
	object_bounds(object_bounds)
{
	evaluate_source_receivers_locations(num_antennas, radius, iteration, location_source, location_receivers);
}

MWInvSystemSetup::MWInvSystemSetup(const MWInvSystemSetup& mw_sys) :
	max_region_length(mw_sys.max_region_length),
	max_region_width(mw_sys.max_region_width),
	sim_region_length(mw_sys.sim_region_length),
	sim_region_width(mw_sys.sim_region_width),
	num_antennas(mw_sys.num_antennas),
	radius(mw_sys.radius),
	iteration(mw_sys.iteration),
	object_type(mw_sys.object_type),
	object_epsilon_r(mw_sys.object_epsilon_r),
	object_sigma(mw_sys.object_sigma),
	object_bounds(mw_sys.object_bounds),
	location_source(mw_sys.location_source),
	location_receivers(mw_sys.location_receivers)
{}

MWInvSystemSetup& MWInvSystemSetup::operator=(const MWInvSystemSetup & mw_sys)
{
	if (this == &mw_sys)
	{
		return *this;
	}

	object_type = mw_sys.object_type;
	object_epsilon_r = mw_sys.object_epsilon_r;
	object_sigma = mw_sys.object_sigma;
	object_bounds = mw_sys.object_bounds;
	location_source = mw_sys.location_source;
	location_receivers = mw_sys.location_receivers;

	return *this;

}

MWInvSystemSetup::~MWInvSystemSetup()
{}

void MWInvSystemSetup::input_system_setup()
{
	double max_region_length;
	double max_region_width;
	double sim_region_length;
	double sim_region_width;
	size_t num_antennas;
	double radius;
	size_t iteration;
	std::string source_type;
	std::vector<double> source_intensities;	
	std::string object_type;
	double object_epsilon_r;
	double object_sigma;
	std::vector<double> object_bounds;

	char response;
	std::cout << "Continue with default system setup? (y/n): ";
	std::cin >> response;
	if (response == 'y')
	{
		return;
	}

	std::cout << "Enter Simulation Region Parameters: " << std::endl;
	std::cout << "Enter the maximum region length: ";
	std::cin >> max_region_length;
	std::cout << "Enter the maximum region width: ";
	std::cin >> max_region_width;
	std::cout << "Enter the simulation region length: ";
	std::cin >> sim_region_length;
	std::cout << "Enter the simulation region width: ";
	std::cin >> sim_region_width;

	std::cout << "Enter Source and Receiver Paramters: " << std::endl;
	std::cout << "Enter the number of total antennas in system: ";
	std::cin >> num_antennas;
	std::cout << "Enter the radius of the antennas: ";
	std::cin >> radius;
	std::cout << "Enter the iteration of the source antenna: ";
	std::cin >> iteration;

	std::cout << "Enter Source Wave Parameters: " << std::endl;
	std::cout << "Enter the source type (sine, square, etc.): ";
	std::cin >> source_type;
	std::cout << "Enter the source intensities: ";
	for (size_t i = 0; i < num_antennas; i++)
	{
		double intensity;
		std::cin >> intensity;
		source_intensities.push_back(intensity);
	}

	std::cout << "Enter the source frequencies: ";
	for (size_t i = 0; i < num_antennas; i++)
	{
		double frequency;
		std::cin >> frequency;
		source_frequencies.push_back(frequency);
	}

	std::cout << "Enter Contrasting Object Parameters: " << std::endl;
	std::cout << "Enter the object type (circle, rectangle, etc.): ";
	std::cin >> object_type;
	std::cout << "Enter the object's relative permittivity: ";
	std::cin >> object_epsilon_r;
	std::cout << "Enter the object's conductivity: ";
	std::cin >> object_sigma;


	if (object_type == "circle")
	{
		double bound;
		std::cout << "Enter the object's radius: ";
		std::cin >> bound;
		object_bounds.push_back(bound);
	}
	else if (object_type == "rectangle")
	{
		double bound;
		std::cout << "Enter the object's length: ";
		std::cin >> bound; object_bounds.push_back(bound);
		std::cout << "Enter the object's width: ";
		std::cin >> bound; object_bounds.push_back(bound);
	}
	else if (object_type == "polygon")
	{
		std::cout << "Enter the number of sides: ";
		size_t num_sides;
		std::cin >> num_sides;
		for (size_t i = 0; i < num_sides; i++)
		{
			std::cout << "Enter the object's side " << i << ": ";
			double side;
			std::cin >> side;
			object_bounds.push_back(side);
		}
	}
	else
	{
		std::cout << "Invalid object type" << std::endl;
	}

	// Calling the overloaded constructor to set the system setup
	MWInvSystemSetup(
		max_region_length, 
		max_region_width, 
		sim_region_length, 
		sim_region_width, 
		num_antennas, 
		radius, 
		iteration,
		source_type,
		source_intensities,
		source_frequencies,
		object_type, 
		object_epsilon_r, 
		object_sigma, 
		object_bounds);
}

// Evaluate the source and receiver locations
void MWInvSystemSetup::evaluate_source_receivers_locations(
	const size_t num_antennas,
	const double radius,
	const size_t iteration,
	Points2D<double>& location_source,
	Points2D<double>& location_receivers)
{
	// Evaluate the source locations
	for (size_t i = 0; i < num_antennas; i++)
	{	
		double angle = 2 * M_PI * i / num_antennas;
		if (i == iteration)
		{
			location_source.x.push_back(radius * cos(angle));
			location_source.y.push_back(radius * sin(angle));
		}
		else
		{
			location_receivers.x.push_back(radius * cos(angle));
			location_receivers.y.push_back(radius * sin(angle));
		}
	}
}

void MWInvSystemSetup::print_system_information()
{
	std::cout << "Simulation Region Parameters: " << std::endl;
	std::cout << "Maximum Region Length: " << max_region_length << std::endl;
	std::cout << "Maximum Region Width: " << max_region_width << std::endl;
	std::cout << "Simulation Region Length: " << sim_region_length << std::endl;
	std::cout << "Simulation Region Width: " << sim_region_width << std::endl;

	std::cout << "Source and Receiver Parameters: " << std::endl;
	std::cout << "Number of Antennas: " << num_antennas << std::endl;
	std::cout << "Radius of Antennas: " << radius << std::endl;
	std::cout << "Iteration of Source Antenna: " << iteration << std::endl;

	std::cout << "Source Wave Parameters: " << std::endl;
	std::cout << "Source Type: " << source_type << std::endl;
	std::cout << "Source Intensities: ";
	for (double i: source_intensities)
	{ std::cout << i << " "; }
	std::cout << std::endl;

	std::cout << "Source Frequencies: ";
	for (double f : source_frequencies)
	{ std::cout << f << " "; }
	std::cout << std::endl;

	std::cout << "Contrasting Object Parameters: " << std::endl;
	std::cout << "Object Type: " << object_type << std::endl;
	std::cout << "Object Relative Permittivity: " << object_epsilon_r << std::endl;
	std::cout << "Object Conductivity: " << object_sigma << std::endl;
	std::cout << "Object Bounds: ";
	for (double b : object_bounds)
	{
		std::cout << b << " ";
	}
	std::cout << std::endl;
}