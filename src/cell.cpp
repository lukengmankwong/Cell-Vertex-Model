#include "cellstructure.h"
#include "tissue.h"
#include "parameters.h"
#include "maths.h"

Cell::Cell(Tissue* T, int id, Halfedge* root)
{
	T_ = T;
	id_ = id;
	root_ = root;
	surface_tension_ = 0;
}

// deletes half edges of cell, sets relevant half edges to have twin nullptrs
void Cell::selfDestroy()
{
	findVerticesAndHalfedges();
	for (int i = 0; i < n_edges_; i++)
	{
		if (halfedges_[i]->twin()->cell() == nullptr) halfedges_[i]->selfDestroy();
		else (halfedges_[i]->removeCell());
	}
	T_->deleteCell(this);
}

void Cell::setRoot(Halfedge* root) { root_ = root; }

void Cell::findVerticesAndHalfedges()
{
	vertices_ = {};
	halfedges_ = {};
	Halfedge* current = root_;
	do {
		vertices_.push_back(current->source());
		halfedges_.push_back(current);
		current = current->next();
	} while (current != root_);
	n_vertices_ = vertices_.size();
	n_edges_ = halfedges_.size();
}

void Cell::boundaryCheck()
{
	Halfedge* current = root_;
	on_boundary_ = false;
	do {
		if (current->twin()->cell() == nullptr) on_boundary_ = true;
		current = current->next();
	} while (current != root_ && !on_boundary_);
}

void Cell::calcArea()
{
	signed_area_ = 0;
	Halfedge* current = root_;
	// shoelace formula
	do {
		double x1 = current->source()->x();
		double y1 = current->source()->y();
		double x2 = current->target()->x();
		double y2 = current->target()->y();
		current = current->next();
		
		signed_area_+= (x1*y2-x2*y1);
		
	} while (current != root_);
	signed_area_ *= 0.5;
	area_ = std::abs(signed_area_);
}

void Cell::calcCentroid()
{
	centroid_x_ = 0;
	centroid_y_ = 0;
	Halfedge* current = root_;
	// shoelace formula taken further to get centroid
	do {
		double x1 = current->source()->x();
		double y1 = current->source()->y();
		double x2 = current->target()->x();
		double y2 = current->target()->y();
		current = current->next();
		
		centroid_x_ += (x1+x2)*(x1*y2-x2*y1);
		centroid_y_ += (y1+y2)*(x1*y2-x2*y1);
		
	} while (current != root_);
	centroid_x_ /= 6*signed_area_;
	centroid_y_ /= 6*signed_area_;
}

void Cell::calcPerimeter()
{
	perimeter_ = 0;
	Halfedge* current = root_;
	do {
		perimeter_ += current->length();
		current = current->next();
	} while (current != root_);
}

void Cell::calcSurfaceTension()
{
	calcArea();
	surface_tension_ = parameter::K_a*(area_-parameter::A_0);
}

void Cell::calcGyration()
{
	G_[0] = 0; G_[1] = 0; G_[2] = 0;
	calcCentroid();
	Halfedge* he_current = root_;
	do {
		Vertex* v = he_current->target();
		double x_v = v->x(); double y_v = v->y();
		G_[0]+=(x_v-centroid_x_)*(x_v-centroid_y_);	
		G_[1]+=(x_v-centroid_x_)*(y_v-centroid_y_);
		G_[2]+=(y_v-centroid_x_)*(y_v-centroid_y_);
		he_current = he_current->next();

	} while(he_current != root_);

	G_[0] /= n_vertices_;
	G_[1] /= n_vertices_;
	G_[2] /= n_vertices_;

	X_ = G_[1];
	Z_ = (G_[0]-G_[2])/2;
	lambda_ = 0.5*( G_[0]+G_[2] + std::sqrt( (G_[0]+G_[2])*(G_[0]+G_[2]) - 4*(G_[0]*G_[2]-G_[1]*G_[1]) ) );
	n_x_ = 1; n_y_ = ((lambda_-G_[0])/G_[1]) / std::sqrt( 1 + ((lambda_-G_[0])/G_[1])*((lambda_-G_[0])/G_[1]));
}

void Cell::calcWinding()
{
	winding_ = 0;
	Halfedge* he_current = root_;
	do {
		if (he_current->cell() && he_current->next()->cell()) winding_ += T_->dTheta(he_current->cell(), he_current->next()->cell());
		he_current = he_current->next();
	} while(he_current != root_);
}

void Cell::output()
{
	std::cout << "id: " << id_ << '\n';
	std::cout << "Vertices: "; for (Vertex* v : vertices_) std::cout << v->id() << " "; std::cout <<'\n';
	std::cout << "n_vertices: " << n_vertices_ << " n_edges: " << n_edges_ << "\n\n";
}

int Cell::id() const { return id_; }
int Cell::n_edges() const { return n_edges_; }
int Cell::n_vertices() const { return n_vertices_; }
const std::vector<Vertex*>& Cell::vertices() const { return vertices_; }
Halfedge* Cell::root() const { return root_; }
bool Cell::on_boundary() const { return on_boundary_; }

double Cell::centroid_x() const { return centroid_x_; }
double Cell::centroid_y() const { return centroid_y_; }
double Cell::area() const { return area_; }
double Cell::signed_area() const { return signed_area_; }
double Cell::perimeter() const { return perimeter_; }
double Cell::surface_tension() const { return surface_tension_; }

double Cell::winding() const { return winding_; }
double Cell::n_x() const { return n_x_; }
double Cell::n_y() const { return n_y_; }
double Cell::X() const { return X_; }
double Cell::Z() const { return Z_; }

