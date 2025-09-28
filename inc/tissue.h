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

	void vertexTranslation();
	void t1Transitions();

public:

	Tissue(const VD& voronoi_diagram, bool (*in)(double, double));
	
	Vertex* createVertex(double x, double y);
	Halfedge* createHalfedge();
	Cell* createCell(Halfedge* root);
	void deleteVertex(Vertex* v);
	void deleteHalfedge(Halfedge* he);
	void deleteCell(Cell* c);

	double dTheta(Cell* c_1, Cell* c_2);
	
	void runTimestep();
	void findDefects();
	void writeCellFile(std::string filename);
	void writeDirectorFile(std::string filename);
	void writeCellDefectFile(std::string filename);

	const std::vector<Vertex*>& vertices() const;
	const std::vector<Halfedge*>& halfedges() const;
	const std::vector<Cell*>& cells() const;

};

#endif // TISSUE_H
