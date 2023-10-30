#include "voronoi.h"

using namespace geovi::algorithm::voronoi_diagram;


void convertPointIndexWithOrigin(geovi::Point2& origin,std::vector<geovi::Point2>& points){
    for(auto& point : points){
        point.x = origin.x - point.x;
        point.y = origin.y - point.y;
    }
}

VoronoiDiagramBuilder::VoronoiDiagramBuilder(NumericalAccuracy numericalAccuracy,Origin originPoint){
    converter = CoordinateSystemConverter(LongitudeBands::band_50);
    converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,originPoint);
    origin = originPoint; 
    multiple = int(numericalAccuracy);
}

void VoronoiDiagramBuilder::build(InputCoordinateSystemType inputType,VoronoiDiagram& voronoi_diagram,Points& points){
   
    if( inputType == CoordinateSystemType::WGS84 ){
        for(auto& point : points){
            converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,point);
        }
    }
    convertPointIndexWithOrigin(origin,points);
    for(auto& point : points){
        point.x = double(int(point.x*multiple));
        point.y = double(int(point.y*multiple));
    }

   // bp::construct_voronoi()

}


