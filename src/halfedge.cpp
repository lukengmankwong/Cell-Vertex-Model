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

void Halfedge::SetSource(Vertex* source)
{
	source_ = source;
}
void Halfedge::SetTarget(Vertex* target)
{
	target_ = target;
}
void Halfedge::SetNext(Halfedge* next)
{
	next_ = next;
}
void Halfedge::SetPrev(Halfedge* prev)
{
	prev_ = prev;
}
void Halfedge::SetCell(Cell* cell)
{
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
void Halfedge::T1Transition()
{
	CalcLength();
	if (length_ > parameter::l_min) return;

	Cell* c_1 = cell_;
	Cell* c_2 = twin_->cell();			// opposite cell
	Cell* c_a = next_->twin()->cell();	// upwards cell
	Cell* c_b = twin_->next()->twin()->cell();	// downwards cell
	if (!(c_1 && c_2 && c_a && c_b) || (c_1->on_boundary() || c_2->on_boundary() || c_a->on_boundary() || c_b->on_boundary())) return;
	
	Vertex* v_a = target_;
	Vertex* v_b = source_;
	Halfedge* he_a_1 = next_->twin();
	Halfedge* he_a_2 = he_a_1->next();
	Halfedge* he_b_1 = twin_->next()->twin();
	Halfedge* he_b_2 = he_b_1->next();
	if (!(v_a && v_b && he_a_1 && he_a_2 && he_b_1 && he_b_2)) return;

	// update vertex positions
	double dx_a = - ( v_b->x() - v_a->x() - (parameter::l_new/length_)*(v_b->y() - v_a->y()) ) / 2;
	double dy_a = - ( v_b->y() - v_a->y() + (parameter::l_new/length_)*(v_b->x() - v_a->x()) ) / 2;
	v_a->Translate(dx_a, dy_a);
	v_b->Translate(-dx_a, -dy_a);

	// update v_1 topology
	v_a->RemoveIncidentEdge(twin_->prev());
	v_a->AddIncidentEdge(prev_);
	// update v_2 topology
	v_b->RemoveIncidentEdge(prev_);
	v_b->AddIncidentEdge(twin_->prev());

	// ensure cell roots are not edges in the transition
	c_1->ChangeRoot(next_);
	c_2->ChangeRoot(twin_->next());

	// update c_1 halfedges
	prev_->SetNext(next_);
	prev_->SetTarget(v_a);
	next_->SetPrev(prev_);
	// update c_2 halfedges
	twin_->prev()->SetNext(twin_->next());
	twin_->prev()->SetTarget(v_b);
	twin_->next()->SetPrev(twin_->prev());
	// update c_a halfedges
	he_a_1->SetNext(twin_);
	he_a_2->SetPrev(twin_);
	he_a_2->SetSource(v_b);
	// update c_b halfedges
	he_b_1->SetNext(this);
	he_b_2->SetPrev(this);
	he_b_2->SetSource(v_a);

	// update twin halfedge
	twin_->SetNext(he_a_2);
	twin_->SetPrev(he_a_1);
	twin_->SetCell(c_a);
	//update this halfedge
	SetNext(he_b_2);
	SetPrev(he_b_1);
	SetCell(c_b);

	c_1->FindVerticesAndEdges();
	c_2->FindVerticesAndEdges();
	c_a->FindVerticesAndEdges();
	c_b->FindVerticesAndEdges();
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

