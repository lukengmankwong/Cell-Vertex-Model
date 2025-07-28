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
void Cell::SelfDestroy()
{
	std::vector<Halfedge*> halfedges;
	Halfedge* current = root_;
	do {
		halfedges.push_back(current);
		current = current->next();
	} while (current != root_);
	
	for (int i = 0; i < halfedges.size(); i++)
	{
		if (halfedges[i]->twin()->cell() == nullptr) halfedges[i]->SelfDestroy();
		else (halfedges[i]->RemoveCell());
	}
	//for (int i = 0; i < halfedges.size(); i++) halfedges[i]->SelfDestroy();
	T_->DeleteCell(this);
}


void Cell::CountEdges()
{
	n_edges_ = 1;
	Halfedge* current = root_->next();
	while (current != root_)
	{
		current = current->next();
		n_edges_++;
	}
}
void Cell::FindVertices()
{
	vertices_ = {root_->target()};

	Halfedge* current = root_->next();
	while (current != root_)
	{
		vertices_.push_back(current->target());
		current = current->next();
	}
	n_vertices_ = vertices_.size();
}

void Cell::BoundaryCheck()
{
	Halfedge* current = root_;
	on_boundary_ = false;
	do {
		if (current->twin()->cell() == nullptr) on_boundary_ = true;
		current = current->next();
	} while (current != root_ && !on_boundary_);
}

// shoelace formula
void Cell::CalcArea()
{
	signed_area_ = 0;
	Halfedge* current = root_;
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
// shoelace formula taken further to get centroid
void Cell::CalcCentroid()
{
	centroid_x_ = 0;
	centroid_y_ = 0;
	Halfedge* current = root_;
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
void Cell::CalcPerimeter()
{
	perimeter_ = 0;
	Halfedge* current = root_;
	do {
		perimeter_ += current->length();
		current = current->next();
	} while (current != root_);
}
void Cell::CalcSurfaceTension()
{
	CalcArea();
	surface_tension_ = parameter::K_a*(area_-parameter::A_0);
}

void Cell::CalcGyration()
{
	G_[0] = 0; G_[1] = 0; G_[2] = 0;
	CalcCentroid();
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
void Cell::CalcWinding()
{
	winding_ = 0;
	Halfedge* he_current = root_;
	do {
		if (he_current->cell() && he_current->next()->cell()) winding_ += T_->DeltaTheta(he_current->cell(), he_current->next()->cell());
		he_current = he_current->next();

	} while(he_current != root_);
}

void Cell::Output()
{
	std::cout << "id: " << id_ << '\n';
	std::cout << "Vertices: "; for (Vertex* v : vertices_) std::cout << v->id() << " "; std::cout <<'\n';
	std::cout << "n_vertices: " << n_vertices_ << " n_edges: " << n_edges_ << "\n\n";
}

int Cell::id() const { return id_; }
int Cell::n_edges() const { return n_edges_; }
int Cell::n_vertices() const { return n_vertices_; }
const std::vector<Vertex*>& Cell::vertices() const { return vertices_; }
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

