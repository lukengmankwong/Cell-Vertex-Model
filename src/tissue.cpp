#include "cellstructure.h"
#include "tissue.h"
#include "diagram.h"
#include "parameters.h"
#include "maths.h"

#include <vector>
#include <unordered_set>
#include <map>
#include <string>
#include <iostream>
#include <fstream>

Tissue::Tissue(const VD& voronoi_diagram, bool (*in)(double, double))
{
	// map CGAL halfedge data structure into custom half edge data structure
	vertex_counter_ = 0;
	halfedge_counter_ = 0;
	cell_counter_ = 0;
	std::map<VD::Vertex_handle, Vertex*> vertex_map;
	std::map<VD::Halfedge_handle, Halfedge*> halfedge_map;
	std::map<VD::Face_handle, Cell*> cell_map;

	// create halfedges
	for (VD::Halfedge_iterator he_it = voronoi_diagram.halfedges_begin(); he_it != voronoi_diagram.halfedges_end(); he_it++)
	{
		if (!he_it->is_unbounded() && !he_it->is_ray())
		{
			Halfedge* he = createHalfedge();
			VD::Halfedge_handle he_h = he_it;
			halfedge_map[he_h] = he;
		}
	}
	// create cells with their root halfedges
	for (VD::Face_iterator f_it = voronoi_diagram.faces_begin(); f_it != voronoi_diagram.faces_end(); f_it++)
	{
		if (!f_it->is_unbounded() && f_it->is_valid())
		{
			VD::Halfedge_handle cell_root = f_it->halfedge();
			Cell* c = createCell(halfedge_map[cell_root]);
			VD::Face_handle f_h = f_it;
			cell_map[f_h] = c;
		}
	}
	// create all vertices
	for (VD::Vertex_iterator v_it = voronoi_diagram.vertices_begin(); v_it != voronoi_diagram.vertices_end(); v_it++)
	{
		if (v_it->is_valid())
		{
			Vertex* v = createVertex(v_it->point().x(), v_it->point().y());
			VD::Vertex_handle v_h = v_it;
			vertex_map[v_h] = v;
		}
	}
	// assign pointers for halfedges
	for (VD::Halfedge_iterator he_it = voronoi_diagram.halfedges_begin(); he_it != voronoi_diagram.halfedges_end(); he_it++)
	{
		if (!he_it->is_unbounded() && !he_it->is_ray())
		{
			Vertex* source = vertex_map[he_it->source()];
			Vertex* target = vertex_map[he_it->target()];
			Halfedge* twin = halfedge_map[he_it->twin()];
			Halfedge* next = halfedge_map[he_it->next()];
			Halfedge* prev = halfedge_map[he_it->previous()];
			Cell* cell = cell_map[he_it->face()];

			VD::Halfedge_handle he_h = he_it;
			halfedge_map[he_h]->setPointers(source, target, twin, next, prev, cell);
		}
	}
	// assign pointers for vertices
	for (VD::Halfedge_iterator he_it = voronoi_diagram.halfedges_begin(); he_it != voronoi_diagram.halfedges_end(); he_it++)
	{
		if (!he_it->is_unbounded() && !he_it->is_ray())
		{
			VD::Vertex_handle v_h = he_it->target();
			Vertex* v = vertex_map[v_h];
			v->addIncidentEdge(halfedge_map[VD::Halfedge_handle(he_it)]);
		}
	}
	// process custom half edge data structure.
	for (Cell* c : cells_) c->findVerticesAndHalfedges();
	for (Halfedge* he : halfedges_) he->calcLength();
	for (Cell* c : cells_)
	{
		c->calcArea();
		c->calcCentroid();
		c->calcPerimeter();
		c->boundaryCheck();
	}
	// remove cells to meet initial condition
	std::vector<Cell*> invalid_cells;
	for (Cell* c : cells_) if (! in(c->centroid_x(), c->centroid_y())) invalid_cells.push_back(c);
	for (Cell* c : invalid_cells) c->selfDestroy();
	for (Cell* c : cells_) c->findVerticesAndHalfedges();
	for (Cell* c : cells_) c->boundaryCheck();

	writeCellFile("cells0.vtk");
	writeDirectorFile("directors0.vtk");
	// sanity check : Euler characteristic should equal 1
	int V = vertices_.size();
	int HE = halfedges_.size(); int E = HE/2;
	int F = cells_.size();
	int Euler = V - E + F;
	std::cout << "V=" << V << "\tE=" << E << "\tF=" << F << '\n';
	std::cout << "V-E+F=" << Euler << '\n';
}

Vertex* Tissue::createVertex(double x, double y)
{
	Vertex* v = new Vertex(this, vertex_counter_, x, y);
	vertices_.push_back(v);
	vertex_counter_++;
	return v;
}

Halfedge* Tissue::createHalfedge()
{
	Halfedge* he = new Halfedge(this, halfedge_counter_);
	halfedges_.push_back(he);
	halfedge_counter_++;
	return he;
}

Cell* Tissue::createCell(Halfedge* root)
{
	Cell* c = new Cell(this, cell_counter_, root);
	cells_.push_back(c);
	cell_counter_++;
	return c;
}

void Tissue::deleteVertex(Vertex* v)
{
	vertices_.erase(std::find(vertices_.begin(), vertices_.end(), v));
	delete v;
}

void Tissue::deleteHalfedge(Halfedge* he)
{
	halfedges_.erase(std::find(halfedges_.begin(), halfedges_.end(), he));
	delete he;
}

void Tissue::deleteCell(Cell* c)
{
	cells_.erase(std::find(cells_.begin(), cells_.end(), c));
	delete c;
}

double Tissue::dTheta(Cell* c_1, Cell* c_2)
{
	double Z_1 = c_1->Z(); double X_1 = c_1->X(); double S_1 = std::sqrt(X_1*X_1+Z_1*Z_1);
	double Z_2 = c_2->Z(); double X_2 = c_2->X(); double S_2 = std::sqrt(X_2*X_2+Z_2*Z_2);
	return std::asin( (Z_1*X_2-X_1*Z_2) / (S_1*S_2) );
}

void Tissue::vertexTranslation()
{
	for (Halfedge* he : halfedges_) he->calcLength();
	for (Cell* c : cells_) 			c->calcPerimeter();
	for (Halfedge* he : halfedges_) he->calcLineTension();

	for (Cell* c : cells_) 			c->calcSurfaceTension();

	for (Vertex* v : vertices_) 	v->updateForce();
	for (Vertex* v : vertices_) 	v->applyForce();
}

void Tissue::t1Transitions()
{
	std::unordered_set<Halfedge*> short_halfedges;
	std::unordered_set<Cell*> cells_of_short_halfedges;
	for (Halfedge* he : halfedges_) he->calcLength();
	for (Halfedge* he : halfedges_) if (he->length() < parameter::l_min && he->cell() && he->twin()->cell() && he->next()->twin()->cell() && he->prev()->twin()->cell())
	{
		short_halfedges.insert(he);
		short_halfedges.erase(he->twin());
	}
	for (Halfedge* he : short_halfedges) he->t1Transition();
}

void Tissue::findDefects()
{
	double tol = 1e-3;

	cell_defects_P_ONE_ = {};
	cell_defects_N_ONE_ = {};
	cell_defects_P_HALF_ = {};
	cell_defects_N_HALF_ = {};

	for (Cell* c : cells_) c->calcGyration();
	for (Cell* c : cells_)
	{
		c->calcWinding();
		double winding = c->winding();
		if (std::fabs(winding - 1) < tol) cell_defects_P_ONE_.push_back(c); 
		else if (std::fabs(winding + 1) < tol) cell_defects_N_ONE_.push_back(c); 
		else if (std::fabs(winding - 0.5) < tol) cell_defects_P_HALF_.push_back(c); 
		else if (std::fabs(winding + 0.5) < tol) cell_defects_N_HALF_.push_back(c); 
	}
}

void Tissue::runTimestep()
{
	vertexTranslation();
	t1Transitions();
	findDefects();
}

void Tissue::writeCellFile(std::string filename)
{
	std::ofstream file(filename);
    file << "# vtk DataFile Version 2.0\nGraph\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << vertices_.size() << " float\n";
	for (Vertex* v : vertices_) { file << v->x() << " " << v->y() << " 0\n"; }

	std::unordered_map<int, int> index_id_map;
    int index = 0;
    for (Vertex* v : vertices_) { index_id_map[v->id()] = index++; }

	int n = cells_.size(); for (Cell* c : cells_) { n += c->n_vertices(); }
	file << "CELLS " << cells_.size() << " " << n << '\n';
	for (Cell* c : cells_)
	{
		file << c->n_vertices() << " ";
		for (Vertex* v : c->vertices()) { file << index_id_map[v->id()] << " "; }
		file << '\n';
	}
	file << "CELL_TYPES " << cells_.size() << '\n';
	for (int c = 0; c < cells_.size(); c++) { file << "7\n"; }

	file.close();
}

void Tissue::writeDirectorFile(std::string filename)
{
	std::ofstream file(filename);
	int n = cells_.size();
	file << "# vtk DataFile Version 2.0\nns\nASCII\nDATASET POLYDATA\nPOINTS " << 2*n << " float\n";
	for (Cell* c : cells_)
	{
		file << (c->centroid_x()-c->n_x()/2) << " " << (c->centroid_y()-c->n_y()/2) << " 0\n";
		file << (c->centroid_x()+c->n_x()/2) << " " << (c->centroid_y()+c->n_y()/2) << " 0\n";
	}
	file << "LINES " << n << " " << 3*n << "\n";
	for (int i = 0; i < n; i++) { file << "2 " << 2*i << " " << 2*i+1 << "\n"; }
	file.close();
}

void Tissue::writeCellDefectFile(std::string filename)
{

}

const std::vector<Vertex*>& Tissue::vertices() const { return vertices_; }
const std::vector<Halfedge*>& Tissue::halfedges() const { return halfedges_; }
const std::vector<Cell*>& Tissue::cells() const { return cells_; }
