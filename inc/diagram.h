#ifndef DIAGRAM_H
#define DIAGRAM_H

#include <vector>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef CGAL::Point_2<K>												  Point;

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Voronoi_diagram_2.h>
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

class Diagram
{
private:

    std::vector<Point> points_;
    VD voronoi_diagram_;

public:

    Diagram(std::vector<Point> points); 

    void createVoronoiDiagram();
    
    const VD& voronoi_diagram() const;

};

#endif // DIAGRAM_H
