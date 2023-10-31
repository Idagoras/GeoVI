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
        double x = double(int(point.x*multiple));
        double y = double(int(point.y*multiple));
        voronoi_diagram.points.push_back(Point2{x,y});
    }
    bp::construct_voronoi(voronoi_diagram.points.begin(),voronoi_diagram.points.end(),&voronoi_diagram.vd);
    

}


VoronoiDiagram::Points VoronoiDiagram::vertices(){
    VoronoiDiagram::Points ps;
    for(bp::voronoi_diagram<double>::const_vertex_iterator it = vd.vertices().begin();
    it != vd.vertices().end();it++){
        ps.push_back(Point2{it->x(),it->y()});
    }
    return ps;
}

