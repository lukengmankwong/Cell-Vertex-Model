#include "diagram.h"

Diagram::Diagram(std::vector<Point> points)
{
    points_ = points;
}

void Diagram::CreateVoronoiDiagram()
{
    DT delauney_tri; 
    delauney_tri.insert(points_.begin(), points_.end());
    voronoi_diagram_ = VD(delauney_tri);
}

const VD& Diagram::VoronoiDiagram() const { return voronoi_diagram_; }
