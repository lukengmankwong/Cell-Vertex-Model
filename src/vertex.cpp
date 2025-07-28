#include "cellstructure.h"
#include "tissue.h"
#include "parameters.h"
#include "maths.h"

Vertex::Vertex(Tissue* T, int id, double x, double y)
{	
	T_ = T;
	id_ = id;
	x_ = x; y_ = y;
	f_x_ = 0; f_y_ = 0;
}

void Vertex::AddIncidentEdge(Halfedge* edge)
{
	incident_edges_.insert(edge);
}
void Vertex::RemoveIncidentEdge(Halfedge* edge)
{
	incident_edges_.erase(edge);
	if (incident_edges_.empty()) SelfDestroy();
}
void Vertex::SelfDestroy()
{
	T_->DeleteVertex(this);
}

void Vertex::UpdateForce()
{
	// line force
	double f_L_x = 0; double f_L_y = 0;
	auto calc_L = [this, &f_L_x, &f_L_y](Vertex* v, Halfedge* he)
	{
		double delta_x = x_ - v->x();
		double delta_y = y_ - v->y();
		double r = d(delta_x, delta_y, 0, 0);
		double dldx = delta_x/r;
		double dldy = delta_y/r;
		f_L_x -= he->line_tension() * dldx;
		f_L_y -= he->line_tension() * dldy;
	};
	for (Halfedge* he : incident_edges_)
	{
		if (he->cell())
		{
			Vertex* v_prev = he->source();
			calc_L(v_prev, he);
			
			Halfedge* he_next = he->next();
			Vertex* v_next = he_next->target();
			calc_L(v_next, he_next);
		}
	}
	
	// area force
	double f_A_x = 0; double f_A_y = 0;
	auto calc_A = [this, &f_A_x, &f_A_y](Halfedge* he)
	{
		if (he->cell())
		{
			Cell* c = he->cell();
			double s = c->signed_area()/c->area();
			double dAdx = 0.5 * (he->next()->target()->x() - he->source()->x()) * s;
			double dAdy = 0.5 * (he->source()->y() - he->next()->target()->y()) * s;
			f_A_x -= c->surface_tension()*dAdx;
			f_A_y -= c->surface_tension()*dAdy;
		}
	};
	for (Halfedge* he : incident_edges_) calc_A(he);

	// update force on vertex
	f_x_ += (f_L_x + f_A_x);
	f_y_ += (f_L_y + f_A_y);
}
void Vertex::ApplyForce()
{
	x_ += parameter::a * f_x_ * parameter::dt;
	y_ += parameter::a * f_y_ * parameter::dt;
}

int Vertex::id() 		const { return id_; }
double Vertex::x() 		const { return x_; }
double Vertex::y() 		const { return y_; }
double Vertex::winding() 	const { return winding_; }

const std::unordered_set<Halfedge*>& Vertex::incident_edges() const { return incident_edges_; }
