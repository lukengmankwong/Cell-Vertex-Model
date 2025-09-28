#include "diagram.h"

Diagram::Diagram(std::vector<Point> points)
{
    points_ = points;
}

void Diagram::createVoronoiDiagram()
{
    DT delauney_tri; 
    delauney_tri.insert(points_.begin(), points_.end());
    voronoi_diagram_ = VD(delauney_tri);
}

const VD& Diagram::voronoi_diagram() const { return voronoi_diagram_; }
