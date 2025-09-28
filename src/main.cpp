#include "cellstructure.h"
#include "tissue.h"
#include "simulation.h"

#include <chrono>

// some initial cell layouts
std::vector<Point> randomPoints(int n) 
{ 
	double A_0 = 1.0;
	std::vector<Point> points;
	double w = A_0*std::sqrt(n);
	for (int i = 0; i < n; i++) points.push_back(Point(w*((static_cast<double>(std::rand())/RAND_MAX)-0.5), w*((static_cast<double>(std::rand())/RAND_MAX)-0.5)));
	return points;
}

// some initial tissue shapes
bool circle(double x, double y) { return x*x + y*y < 1000; }
bool annulus(double x, double y) { return x*x + y*y < 1000 && x*x + y*y > 100; }

int main()
{
	int n_initial_sites = 10000;
	int n_timesteps = 1000;
	int output_interval = 1;

	// initial condition
	std::cout << "CREATING VORONOI DIAGRAM\n";
	std::vector<Point> points = randomPoints(n_initial_sites);
	Diagram diagram = Diagram(points);
	diagram.createVoronoiDiagram();

	std::cout << "INITIALISING TISSUE\n";
	auto init_start = std::chrono::high_resolution_clock::now();
	Tissue T(diagram.voronoi_diagram(), circle);
	auto init_end = std::chrono::high_resolution_clock::now();
	auto init_duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(init_end - init_start);
	std::cout << "INITIALISED IN " << init_duration_ms.count() << " ms\n";

	// simulation
	Simulation simulation(&T);
	std::cout << "RUNNING SIMULATION\n";
    auto sim_start = std::chrono::high_resolution_clock::now();
	simulation.run(n_timesteps, output_interval);
	auto sim_end = std::chrono::high_resolution_clock::now();
    auto sim_duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(sim_end - sim_start);
	std::cout << "SIMULATION RAN IN " << sim_duration_ms.count() << " ms\n";
	
	return 1;
}

// To do:
// cell extrusion
// cell division
