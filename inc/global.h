#ifndef GLOBALS_H
#define GLOBALS_H

#include <utility>
#include <unordered_map>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef CGAL::Vector_2<K>                                                   Vec;
typedef CGAL::Point_2<K>												  Point;

#include "vertex.h"
#include "edge.h"
#include "cell.h"

class Vertex;
class Edge;
class Cell;

//constants

extern const double pi;

//parameters
extern const int cell_count;

extern const double dt;
extern const int timesteps;

extern const double A_0; //characteristic cell area

extern const double k_A; //spring constant for surface tension
extern const double k_L; //spring constant for edge tension
extern const double T_l_0; //independent edge tension

extern const double a;

class Global //singleton
{
private:
	
	Global();
	~Global();

	std::unordered_map<int, Vertex> v_map; int v_c_;
	std::unordered_map<int, Edge> e_map; int e_c_;
	std::unordered_map<int, Cell> c_map; int c_c_;
	
	int timestep;
	std::vector<std::vector<int>> cell_defects;
	std::vector<std::vector<int>> vertex_defects;
	
	void extrusion();
	void division();

	
public:
	void T1();
	Global(const Global&) = delete;
    Global& operator=(const Global&) = delete;
    static Global& get() 
    {
        static Global instance;
        return instance;
    }
    
    void nextStep();
    const int Step() const;
    void addDefects();
    const std::vector<int>& cellStepDefects() const;
    const std::vector<int>& vertexStepDefects() const;
    
    std::unordered_map<int, Vertex>& vertexMap();
    std::unordered_map<int, Edge>& edgeMap();
    std::unordered_map<int, Cell>& cellMap();
    
    Vertex& vert(int v);
    Edge& edge(int e);
    Cell& cell(int c);
    
    const int v_c() const; 
    const int e_c() const;
    const int c_c() const;

	const int createVertex(Point r);
	const int createEdge(int v1, int v2);
	const int createCell(std::vector<int>& vertices, std::vector<int>& edges);
	
	void cellNewVertex(int c, int v, int i); //add new vertex to cell vertices at index i
	void cellNewEdge(int c, int e, int i); //add new edge to cell edges at index i
	
	void cellRemoveVertex(int c, int v);
	int cellRemoveEdge(int c, int e);
	
	void cellExchangeVertex(int c, int v_old, int v_new);
	
	void destroyVertex(int v);
	void destroyEdge(int e);
	void destroyCell(int c);
	
	const bool commonEdge(int c1, int c2) const;
	
	void transitions();
	void run();
	
	const double D_angle(int c_i, int c_j) const; 
};

#endif // GLOBALS_H
