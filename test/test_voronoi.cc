#include "GeoVI/voronoi.h"
#include <iostream>
#include <vector>
using namespace geovi::algorithm::voronoi_diagram;

int main(){
    VoronoiDiagramBuilder builder(geovi::NumericalAccuracy::meter,geovi::Point2{0,0});
    VoronoiDiagram vd;
    std::vector<geovi::Point2> points;
    points.push_back(geovi::Point2{1,0});
    points.push_back(geovi::Point2{1,6});
    points.push_back(geovi::Point2{10,4});
    points.push_back(geovi::Point2{7,5});
    points.push_back(geovi::Point2{2,9});
    builder.build(VoronoiDiagramBuilder::InputCoordinateSystemType::UTM,vd,points);
    VoronoiDiagram::Points ps  = vd.vertices();
    for (auto p : ps){
        std::cout << "vertex coordinate x = " << p.x << " vertex coordinate y " << p.y << std::endl;
    }
    
    return 0;
}