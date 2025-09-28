#ifndef TISSUE_H
#define TISSUE_H

#include "diagram.h"

#include <string>
#include <vector>

class Vertex;
class Halfedge;
class Cell;

class Tissue
{
private:

	int vertex_counter_;
	int halfedge_counter_;
	int cell_counter_;

	std::vector<Vertex*> vertices_;
	std::vector<Halfedge*> halfedges_;
	std::vector<Cell*> cells_;

	std::vector<Cell*> cell_defects_P_ONE_;
	std::vector<Cell*> cell_defects_N_ONE_;
	std::vector<Cell*> cell_defects_P_HALF_;
	std::vector<Cell*> cell_defects_N_HALF_;

	void VertexTranslation();
	void T1Transitions();

public:

	Tissue(const VD& voronoi_diagram, bool (*in)(double, double));
	
	Vertex* CreateVertex(double x, double y);
	Halfedge* CreateHalfedge();
	Cell* CreateCell(Halfedge* root);
	void DeleteVertex(Vertex* v);
	void DeleteHalfedge(Halfedge* he);
	void DeleteCell(Cell* c);

	double DeltaTheta(Cell* c_1, Cell* c_2);
	
	void RunTimestep();
	void FindDefects();
	void WriteCellFile(std::string filename);
	void WriteDirectorFile(std::string filename);
	void WriteCellDefectFile(std::string filename);

	const std::vector<Vertex*>& vertices() const;
	const std::vector<Halfedge*>& halfedges() const;
	const std::vector<Cell*>& cells() const;

};

#endif // TISSUE_H
