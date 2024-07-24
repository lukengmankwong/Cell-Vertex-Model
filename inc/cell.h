#ifndef CELL_H
#define CELL_H

#include <cmath>
#include <array>
#include <vector>
#include <unordered_map>

#include <CGAL/Origin.h>

#include "vertex.h"
#include "edge.h"

class Cell
{
private:

	const int id;
	
    std::unordered_set<int> vertex_keys;
    std::unordered_set<int> edge_keys;

    Point centroid;
    double area;
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double lambda; //gyration tensor largest eigenvalue
    Vec director;
    
    std::unordered_map<int, Vertex>* vertex_map_ptr; int* vertex_counter_ptr;
    std::unordered_map<int, Edge>* edge_map_ptr; int* edge_counter_ptr;
    std::unordered_map<int, Cell>* cell_map_ptr; int* cell_counter_ptr;
    
public:

    Cell(int* vertex_counter_ptr, std::unordered_map<int, Vertex>* vertex_map_ptr, std::vector<int>& vertex_keys,
     int* edge_counter_ptr, std::unordered_map<int, Edge>* edge_map_ptr, std::vector<int>& edge_keys, 
     int* cell_counter_ptr, std::unordered_map<int, Cell>* cell_map_ptr);
     
    ~Cell();
      
	bool removeEdge(int edge_id);
	bool removeVertex(int vertex_id);
        
    void removeEdges();
    void removeVertices();
    void extrude();

    void calcCentroid();
    void calcArea();
    void calcG();
    void calcLambda();

    const int getID() const;
    const Point& getCentroid() const;
    double getArea() const;
    const std::unordered_set<int>& getVertices() const;
    const std::unordered_set<int>& getEdges() const;

};

#endif // CELL_H
