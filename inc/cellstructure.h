#ifndef CELLSTRUCTURE_H
#define CELLSTRUCTURE_H

#include <vector>
#include <unordered_set>

class Vertex;
class Halfedge;
class Cell;
class Tissue;
class Gyration;

class Vertex
{
private:
	
	int id_;
	double x_, y_;
	double f_x_, f_y_;
	double winding_;

	Tissue* T_;
	std::unordered_set<Halfedge*> incident_edges_;
	

public:

	Vertex(Tissue* T, int id, double x, double y);

	void addIncidentEdge(Halfedge* edge);
	void removeIncidentEdge(Halfedge* edge);
	void selfDestroy();

	void translate(double dx, double dy);
	void updateForce();
    void applyForce();

	int id() const;
	double x() const; double y() const;
	double winding() const;

	const std::unordered_set<Halfedge*>& incident_edges() const;

};

class Halfedge
{
private:

	int id_;
	double length_;
	double line_tension_;

	Tissue* T_;
	Vertex* source_;
	Vertex* target_;
	Halfedge* twin_;
	Halfedge* next_;
	Halfedge* prev_;
	Cell* cell_;


public:

	Halfedge(Tissue* T, int id);

	void setPointers(Vertex* source, Vertex* target, Halfedge* twin, Halfedge* next, Halfedge* prev, Cell* cell);
	void setSource(Vertex* source);
	void setTarget(Vertex* target);
	void setNext(Halfedge* next);
	void setPrev(Halfedge* prev);
	void setCell(Cell* cell);
	void removeCell();
	void twinDestroy();
	void selfDestroy();
	void t1Transition();

	void calcLength();
	void calcLineTension();

	int id() const;
	double length() const;
	double line_tension() const;

	Vertex* source() const;
	Vertex* target() const;
	Halfedge* twin() const;
	Halfedge* next() const;
	Halfedge* prev() const;
	Cell* cell() const;

};

class Cell
{
private:

	int id_;
	int n_edges_;
	int n_vertices_;
	bool on_boundary_;

	double signed_area_;
	double area_;
	double centroid_x_, centroid_y_;
	double perimeter_;
	double surface_tension_;

	double G_[3]; // gyration tensor
	double lambda_;
	double X_; double Z_;
	double n_x_; double n_y_;
	double winding_;

	Tissue* T_;
	Halfedge* root_;
	std::vector<Vertex*> vertices_;
	std::vector<Halfedge*> halfedges_;

public:

	Cell(Tissue* T, int id, Halfedge* root);

	void selfDestroy();
	void setRoot(Halfedge* root);

	void findVerticesAndHalfedges();
	void boundaryCheck();
	
	void calcArea();
	void calcCentroid();
	void calcPerimeter();
	void calcSurfaceTension();

	void calcGyration();
	void calcWinding();

	int id() const;	
	int n_edges() const;
	int n_vertices() const;
	const std::vector<Vertex*>& vertices() const;
	Halfedge* root() const;
	bool on_boundary() const;

	double area() const;
	double signed_area() const;
	double centroid_x() const;
	double centroid_y() const;
	double perimeter() const;
	double surface_tension() const;

	double winding() const;
	double n_x() const;
	double n_y() const;
	double X() const;
	double Z() const;
	
	void output();

};

#endif // CELLSTRUCTURE_H
