#include "voronoi.h"
#include <queue>

using namespace geovi::algorithm::voronoi_diagram;


void convertPointIndexWithOrigin(geovi::Point2& origin,std::vector<geovi::Point2>& points){
    for(auto& point : points){
        point.x = origin.x - point.x;
        point.y = origin.y - point.y;
    }
}

void BFS(VoronoiDiagram::Matrix& matrix,VoronoiDiagram::VD& vd){
    std::queue<VoronoiDiagram::VD::cell_type*> cqueue;
    VoronoiDiagram::Matrix visit(matrix.size1(),matrix.size2());
    for(auto it = vd.cells().begin();it!=vd.cells().end();it++){
        visit.insert_element(it->source_index(),it->source_index(),1);
        std::queue<VoronoiDiagram::VD::cell_type*> empty;
        std::swap(empty,cqueue);
        auto source = *it;
        cqueue.push(&source);
        while( !cqueue.empty()){
            auto cell = cqueue.front();
            if (!cell->is_degenerate() ){
                auto edge = cell->incident_edge();
                do{
                    auto neighbor = edge->twin()->cell();
                    if(!visit.at_element(it->source_index(),neighbor->source_index())){
                        matrix.insert_element(it->source_index(),neighbor->source_index(),matrix.at_element(it->source_index(),cell->source_index())+1);
                        cqueue.push(neighbor);
                        visit.insert_element(it->source_index(),neighbor->source_index(),1);
                    }
                } while( edge != cell->incident_edge());
            }
            cqueue.pop();
        }
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

VoronoiDiagram::Matrix VoronoiDiagram::shortestPathBetweenCells(){
    VoronoiDiagram::Matrix m(vd.num_cells(),vd.num_cells());
    BFS(m,vd);
    return m;
}

VoronoiDiagram::Points VoronoiDiagram::sites(){
    return points;
}

VoronoiDiagram::Segments VoronoiDiagram::finiteEdges(){
    VoronoiDiagram::Segments segments;
    for(auto it = vd.edges().begin(); it != vd.edges().end(); it ++ ){
        if( it->is_finite()){
            Segment segment = Segment(it->vertex0()->x(),it->vertex0()->y(),it->vertex1()->x(),it->vertex1()->y());
            segments.push_back(segment);
        }  
    }
    return segments;
}

VoronoiDiagram::Lines VoronoiDiagram::infiniteEdges(){
    VoronoiDiagram::Lines lines;
    for(auto it = vd.edges().begin(); it != vd.edges().end() ; it ++ ){
        if ( it -> is_infinite()){
            Line line = Line();
            if ( it->vertex0() == NULL ){
                line.p0 = Point2(it->vertex1()->x(),it->vertex1()->y());
                auto cell = it->cell();
                auto neighbor = it->twin()->cell();
                line.slope = (points[cell->source_index()].x-points[neighbor->source_index()].x)/(points[cell->source_index()].y-points[neighbor->source_index()].y);
                line.dir = Direction::in;
            }
            if ( it->vertex1() == NULL ){
                line.p0 = Point2(it->vertex0()->x(),it->vertex0()->y());
                auto cell = it->cell();
                auto neighbor = it->twin()->cell();
                line.slope = (points[cell->source_index()].x-points[neighbor->source_index()].x)/(points[cell->source_index()].y-points[neighbor->source_index()].y);
                line.dir = Direction::out;
            }
            lines.push_back(line);
        }
    }
    return lines;
}