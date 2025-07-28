#include "cellstructure.h"
#include "tissue.h"
#include "parameters.h"
#include "maths.h"

Halfedge::Halfedge(Tissue* T, int id)//, int origin_id, int twin_id, int next_id, int prev_id, int cell_id)
{
	T_ = T;
	id_ = id;

	source_ = nullptr;
	twin_ = nullptr;
	next_ = nullptr;
	prev_ = nullptr;
	cell_ = nullptr;

	line_tension_ = 0;
	
}

void Halfedge::SetPointers(Vertex* source, Vertex* target, Halfedge* twin, Halfedge* next, Halfedge* prev, Cell* cell)
{
	source_ = source;
	target_ = target;

	twin_ = twin;
	next_ = next;
	prev_ = prev;
	
	cell_ = cell;
}

void Halfedge::RemoveCell()
{
	cell_ = nullptr;
}

void Halfedge::TwinDestroy()
{
	twin_->target()->RemoveIncidentEdge(twin_);
	T_->DeleteHalfedge(twin_);
}
void Halfedge::SelfDestroy()
{
	TwinDestroy();
	target_->RemoveIncidentEdge(this);
	T_->DeleteHalfedge(this);
}

void Halfedge::CalcLength()
{
	double x1 = source_->x();
	double y1 = source_->y();
	double x2 = target_->x();
	double y2 = target_->y();

	length_ = d(x1,y1, x2,y2);
}
void Halfedge::CalcLineTension()
{
	if (cell_ && twin_->cell()) line_tension_ = parameter::Lambda + parameter::Gamma*(cell_->perimeter() + twin_->cell()->perimeter());
}

int Halfedge::id() const { return id_; }
double Halfedge::length() const { return length_; }
double Halfedge::line_tension() const { return line_tension_; }

Vertex* Halfedge::source() const { return source_; }
Vertex* Halfedge::target() const { return target_; }
Halfedge* Halfedge::twin() const { return twin_; }
Halfedge* Halfedge::next() const { return next_; }
Halfedge* Halfedge::prev() const { return prev_; }
Cell* Halfedge::cell() const { return cell_; }

