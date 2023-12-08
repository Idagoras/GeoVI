#include "GeoVI/voronoi.h"
#include <iostream>
#include <vector>
#include <limits>
#include <ctime>
#include <random>



int main(){
    using namespace geovi::algorithm::voronoi_diagram;
    using namespace geovi;

    double d1 = DistanceCalculator::euclidDistance2D({7.7,10},{4,8});
    double d2 = DistanceCalculator::euclidDistance2D({7.7,10},{9,6});
    if(DistanceCalculator::double_distance_equal(d1,d2)){
        std::cout << "is equal" << std::endl;
    }else{
        std::cout << "not equal "<< std::abs(d1-d2) << std::endl;
    }

    std::vector<Point2> rect_boundary;
    std::vector<Point2> sites;
    rect_boundary.emplace_back(0,0);
    rect_boundary.emplace_back(0,10);
    rect_boundary.emplace_back(10,10);
    rect_boundary.emplace_back(10,0);

    std::default_random_engine e;
    std::uniform_int_distribution<int> u(0,5);
    e.seed(std::time(0));

    for(int i = 0 ; i< 5 ; i++)
        sites.emplace_back(u(e),u(e));
        
   // sites.emplace_back(3,4);
   // sites.emplace_back(6,2);
    //sites.emplace_back(6,6);
   // sites.emplace_back(4,8);
   // sites.emplace_back(9,6);
    for(auto s : sites){
        std::cout << "s.x = " << s.x << " y = " << s.y <<std::endl;
    }


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