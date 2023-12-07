#include "GeoVI/voronoi.h"
#include <iostream>
#include <vector>
#include <limits>


int main(){
    using namespace geovi::algorithm::voronoi_diagram;
    using namespace geovi;
    std::vector<Point2> rect_boundary;
    std::vector<Point2> sites;
    rect_boundary.emplace_back(0,0);
    rect_boundary.emplace_back(0,10);
    rect_boundary.emplace_back(10,10);
    rect_boundary.emplace_back(10,0);

    sites.emplace_back(3,4);
    sites.emplace_back(6,2);
    sites.emplace_back(6,6);
    sites.emplace_back(4,8);
    sites.emplace_back(9,6);



    VoronoiDiagramBuilder builder(NumericalAccuracy::meter,rect_boundary);
    VoronoiDiagram voronoi_diagram;
    builder.build(voronoi_diagram,sites);

   for(int i = 0; i < sites.size() ; ++ i ){
       std::cout << "cell index " << i << " site x = " << sites[i].x << " y = " << sites[i].y << std::endl;
       auto cell_polygon = voronoi_diagram.cellPolygon(i);

       std::cout << "polygon:" << std::endl;
       for(auto p : cell_polygon){
           std::cout << "vertex p.x = " << p.x << " p.y = " << p.y <<std::endl;
       }
       std::cout << std::endl;
   }
    
    return 0;
}