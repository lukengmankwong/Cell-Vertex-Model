#ifndef VERTEX_H
#define VERTEX_H

#include <utility>
#include <unordered_set>

#include "global.h"
class Global;

class Vertex
{
private:

	Global* g;
    const int id;
    Point r; Vec dr;
    Vec force;
    
    std::unordered_set<int> cell_contacts;
    std::unordered_set<int> edge_contacts;
    
    Vec calcSurfaceForce();
    Vec calcLineForce();

public:
    Vertex(Global* g, int id, Point r);
    bool operator==(const Vertex& other) const;
    
    const Point& getR() const;
    const int getID() const;
    const std::unordered_set<int>& getCellContacts() const;

    void addCellContact(int cell_id);
    void removeCellContact(int cell_id);
  
    void addEdgeContact(int edge_id);
    void removeEdgeContact(int edge_id);

    void calcForce();
    void addNoise();
    void applyForce();

	void destroy();
	
};

#endif // VERTEX_H
