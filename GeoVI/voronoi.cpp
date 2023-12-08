#include "voronoi.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <deque>
#include <limits>





using namespace geovi::algorithm::voronoi_diagram;
using namespace geovi;
using namespace boost::geometry;




void convertPointIndexWithOrigin(geovi::Point2& origin,std::vector<geovi::Point2>& points){
    for(auto& point : points){
        point.x = origin.x - point.x;
        point.y = origin.y - point.y;
    }
}

void BFS(VoronoiDiagram::Matrix& matrix,const VoronoiDiagram::VD& vd){
    std::queue<VoronoiDiagram::VD::cell_type*> cell_queue;
    boost::numeric::ublas::matrix<bool> visit(matrix.size1(),matrix.size2());
    for(auto it = vd.cells().begin();it!=vd.cells().end();it++){
        visit.insert_element(it->source_index(),it->source_index(),true);
        std::queue<VoronoiDiagram::VD::cell_type*> empty;
        std::swap(empty, cell_queue);
        auto source = *it;
        cell_queue.push(&source);
        while( !cell_queue.empty()){
            auto cell = cell_queue.front();
            if (!cell->is_degenerate() ){
                auto edge = cell->incident_edge();
                do{
                    auto neighbor = edge->twin()->cell();
                    if(!visit.at_element(it->source_index(),neighbor->source_index())){
                        matrix.insert_element(it->source_index(),neighbor->source_index(),matrix.at_element(it->source_index(),cell->source_index())+1);
                        cell_queue.push(neighbor);
                        visit.insert_element(it->source_index(),neighbor->source_index(),true);
                    }
                } while( edge != cell->incident_edge());
            }
            cell_queue.pop();
        }
    } 
}

VoronoiDiagramBuilder::VoronoiDiagramBuilder(NumericalAccuracy numericalAccuracy,const std::vector<Point2>& rect_boundary,Origin originPoint):
m_multiple(int(numericalAccuracy)),
m_origin(originPoint),
m_min_x(std::numeric_limits<double>::max()),
m_min_y(std::numeric_limits<double>::max()),
m_max_x(std::numeric_limits<double>::min()),
m_max_y(std::numeric_limits<double>::min())
{
    for(auto vertex : rect_boundary){
        m_min_x = vertex.x < m_min_x ? vertex.x : m_min_x;
        m_min_y = vertex.y < m_min_y ? vertex.y : m_min_y;
        m_max_x = vertex.x > m_max_x ? vertex.x : m_max_x;
        m_max_y = vertex.y > m_max_y ? vertex.y : m_max_y;
    }
}

void VoronoiDiagramBuilder::build(VoronoiDiagram& voronoi_diagram,Points& points){
    for(auto& point : points){
        auto x = double(int(point.x * m_multiple));
        auto y = double(int(point.y * m_multiple));
        voronoi_diagram.m_points.push_back(Point2{x, y});
    }
    bp::construct_voronoi(voronoi_diagram.m_points.begin(), voronoi_diagram.m_points.end(), &voronoi_diagram.m_vd);
    std::vector<Point2> boundary;
    boundary.emplace_back(m_min_x,m_min_y);
    boundary.emplace_back(m_min_x,m_max_y);
    boundary.emplace_back(m_max_x,m_max_y);
    boundary.emplace_back(m_max_x,m_min_y);
    boundary.emplace_back(m_min_x,m_min_y);

    voronoi_diagram.m_boundary = boundary;

}


VoronoiDiagram::Points VoronoiDiagram::vertices() const{
    VoronoiDiagram::Points ps;
    for(auto it = m_vd.vertices().begin();
        it != m_vd.vertices().end(); it++){
        ps.push_back(Point2{it->x(),it->y()});
    }
    return ps;
}

VoronoiDiagram::Matrix VoronoiDiagram::shortestPathBetweenCells() const{
    VoronoiDiagram::Matrix m(m_vd.num_cells(), m_vd.num_cells());
    BFS(m, m_vd);
    return m;
}

VoronoiDiagram::Points VoronoiDiagram::sites() const{
    return m_points;
}

VoronoiDiagram::Segments VoronoiDiagram::finiteEdges() const{
    VoronoiDiagram::Segments segments;
    for(auto it = m_vd.edges().begin(); it != m_vd.edges().end(); it ++ ){
        if( it->is_finite()){
            Segment segment = Segment(it->vertex0()->x(),it->vertex0()->y(),it->vertex1()->x(),it->vertex1()->y());
            segments.push_back(segment);
        }  
    }
    return segments;
}

VoronoiDiagram::Lines VoronoiDiagram::infiniteEdges() const{
    VoronoiDiagram::Lines lines;
    for(auto it = m_vd.edges().begin(); it != m_vd.edges().end() ; it ++ ){
        if ( it -> is_infinite()){
            Line line = Line();
            if ( it->vertex0() == nullptr ){
                line.p0 = Point2(it->vertex1()->x(),it->vertex1()->y());
                auto cell = it->cell();
                auto neighbor = it->twin()->cell();
                line.slope = -(m_points[cell->source_index()].x - m_points[neighbor->source_index()].x) / (m_points[cell->source_index()].y - m_points[neighbor->source_index()].y);
                line.dir = Direction::in;
            }
            if ( it->vertex1() == nullptr ){
                line.p0 = Point2(it->vertex0()->x(),it->vertex0()->y());
                auto cell = it->cell();
                auto neighbor = it->twin()->cell();
                line.slope = -(m_points[cell->source_index()].x - m_points[neighbor->source_index()].x) / (m_points[cell->source_index()].y - m_points[neighbor->source_index()].y);
                line.dir = Direction::out;
            }
            lines.push_back(line);
        }
    }
    return lines;
}

std::pair<Point2,VoronoiDiagram::CellIndex> VoronoiDiagram::cellIncludePoint(Point2 p) const{
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    double max = std::numeric_limits<double>::max();
    Point2 nearest_cell;
    CellIndex nearest_cell_index = -1;
    CellIndex index = 0 ;
    for(auto point : m_points){
        double distance = boost::geometry::distance(point_xy{p.x,p.y},point_xy{point.x,point.y});
        if( distance < max ){
            max = distance;
            nearest_cell = point;
            nearest_cell_index = index;
        }
        index ++;
    }
    return std::make_pair(nearest_cell,nearest_cell_index);
}

std::vector<VoronoiDiagram::CellIndex> VoronoiDiagram::neighbors(VoronoiDiagram::CellIndex index) const{
    std::vector<VoronoiDiagram::CellIndex> neighbors_indexes;
    for(auto it = m_vd.cells().begin(); it != m_vd.cells().end(); it ++ ){
        if( it -> source_index() == index ){
            auto cell = *it;
            auto edge = cell.incident_edge();
            do{
                VoronoiDiagram::CellIndex neighbor_index = edge->twin()->cell()->source_index();
                neighbors_indexes.push_back(neighbor_index);
                edge = edge->next();
            }while(edge != cell.incident_edge());
        }
    }
    return neighbors_indexes;
}

void intersection_between_line_and_rectangle(const Point2& p,double slope,std::vector<Point2>& intersection_points,const std::vector<Point2>& boundary){
    double x_min = std::numeric_limits<double>::max(),x_max = std::numeric_limits<double>::min(),
            y_min = std::numeric_limits<double>::max() ,y_max = std::numeric_limits<double>::min();
    for(auto point : boundary){
        x_min = point.x < x_min ? point.x : x_min;
        y_min = point.y < y_min ? point.y : y_min;
        x_max = point.x > x_max ? point.x : x_max;
        y_max = point.y > y_max ? point.y : y_max;
    }
    if(slope == std::numeric_limits<double>::infinity()){
        intersection_points.emplace_back(p.x,y_max);
        intersection_points.emplace_back(p.x,y_min);

    }else{
        Point2 intersection_left_x,intersection_right_x,intersection_top_y,intersection_bottom_y,intersection_left,intersection_right;
        if(slope == 0){
            intersection_left = {x_min,p.y};
            intersection_right = {x_max,p.y};
        }else{
            intersection_left_x = {x_min,slope * (x_min - p.x) + p.y};
            intersection_right_x = {x_max,slope * (x_max - p.x) + p.y};
            intersection_top_y = {(y_max-p.y)/slope + p.x,y_max};
            intersection_bottom_y = {(y_min-p.y)/slope + p.x,y_min};
            if(intersection_left_x.y > y_max || intersection_left_x.y < y_min)
                intersection_left = intersection_top_y.x > intersection_bottom_y.x ? intersection_bottom_y : intersection_top_y;
            else
                intersection_left = intersection_left_x;

            if(intersection_right_x.y > y_max || intersection_right_x.y < y_min)
                intersection_right = intersection_top_y.x > intersection_bottom_y.x ? intersection_top_y : intersection_bottom_y;
            else
                intersection_right = intersection_right_x;
        }

        intersection_points.push_back(intersection_left);
        intersection_points.push_back(intersection_right);
    }

}

std::vector<Point2> generate_polygon_with_rect_vertexes_and_line(const std::vector<Point2>& rect_vertexes,const std::vector<Point2>& line){
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::multi_point<point_xy> multi_points;
    typedef boost::geometry::model::polygon<point_xy> polygon;
    multi_points points;
    for(auto& v : rect_vertexes){
        points.emplace_back(v.x,v.y);
    }
    for(auto& v : line){
        points.emplace_back(v.x,v.y);
    }
    polygon hull;
    boost::geometry::convex_hull(points,hull);
    auto outers = hull.outer();
    std::vector<Point2> results;
    for(auto& v : outers){
        results.emplace_back(v.x(),v.y());
    }
    return results;
}


void intersection_polygon(const std::vector<Point2>& rect,const std::vector<Point2>& polygon,std::vector<Point2>& out){
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::polygon<point_xy> b_polygon;
    b_polygon b_rect,b_po;
    for(auto& v : rect){
        boost::geometry::append(b_rect,point_xy{v.x,v.y});
    }

    for(auto& v : polygon){
        boost::geometry::append(b_po,point_xy{v.x,v.y});
    }

    std::deque<b_polygon> b_out;
    boost::geometry::intersection(b_rect,b_po,b_out);
    std::cout << boost::geometry::wkt(b_out[0]) << std::endl;
    out.clear();
    for(auto v : b_out[0].outer()){
        out.emplace_back(v.x(),v.y());
    }
}


std::vector<VoronoiDiagram::CellIndex> VoronoiDiagram::cells_include_or_on_its_edge(Point2& p) const{
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    std::vector<CellIndex> results;
    double min_distance = std::numeric_limits<double>::max();
    uint64_t index = 0;
    for(auto site : m_points){
        double distance = DistanceCalculator::euclidDistance2D(p,site);
        std::cout << "site " << index <<" to point distance = " << distance <<std::endl;
        if(DistanceCalculator::double_distance_equal(min_distance,distance) > 0){
            results.clear();
            results.push_back(index);
            min_distance = distance;
        }
        else if( DistanceCalculator::double_distance_equal(min_distance,distance) == 0)
            results.push_back(index);
        ++ index;
    }
    for( auto r : results){
        std::cout << "include index = " << r << std::endl;
    }
    return results;
}

boost::polygon::voronoi_cell<double>* get_cell_obj(uint64_t index,const boost::polygon::voronoi_diagram<double>& vd){
    auto cells = vd.cells();
    for(auto& cell : cells){
        if( cell.source_index() == index){
            return &cell;
        }
    }
    return nullptr;
}

std::vector<Point2> VoronoiDiagram::cellPolygon(VoronoiDiagram::CellIndex index) const{

    std::vector<Point2> polygon_points;
    auto cells = m_vd.cells();
    boost::polygon::voronoi_cell<double>* cell;
    for(auto& c : cells){
        if( c.source_index() == index){
            cell = &c;
        }
    }
    auto cell_point = m_points[index];
    auto begin = cell->incident_edge();
    auto e = begin;

    // 求位于生成点单元内的矩形边界顶点
    std::vector<Point2> boundary_vertexes_in_cell;
    for(auto b_vertex : m_boundary){
        auto results = cells_include_or_on_its_edge(b_vertex);

        auto find_result = std::find(results.begin(),results.end(),index);
        if(find_result != results.end()){
            boundary_vertexes_in_cell.push_back(b_vertex);
        }

    }
    do{
        if(e->is_finite()){

            polygon_points.emplace_back(e->vertex0()->x(),e->vertex0()->y());
            polygon_points.emplace_back(e->vertex1()->x(),e->vertex1()->y());
        }else{
            // 求存在的一个端点
            auto v = e->vertex0() == nullptr ? e->vertex1() : e->vertex0();
            // 求相邻的生成点指针
            auto neighbor = e->twin()->cell();
            // 求相邻的生成点坐标
            auto neighbor_point = m_points[neighbor->source_index()];
            // 求相邻生成点坐标和生成点坐标的中点，该点位于这条无限边上
            Point2 mid_point = {(cell_point.x+neighbor_point.x)/2,(cell_point.y+neighbor_point.y)/2};
            // 求该条无限边的斜率
            double slope = cell_point.y - neighbor_point.y == 0 ? std::numeric_limits<double>::infinity() : -(cell_point.x-neighbor_point.x)/(cell_point.y-neighbor_point.y);
            // 求无线长的边和边界框的两个交点，返回值的首个元素相比第二个元素位置更左，也有可能有相等的x坐标
            std::vector<Point2> intersection_points ;
            intersection_between_line_and_rectangle(mid_point,slope,intersection_points,m_boundary);
            std::vector<Point2> intersection_point_in_cell;
            for(auto inter_p : intersection_points){
                std::cout << "intersection" << std::endl<<std::endl; 
                auto results = cells_include_or_on_its_edge(inter_p);

                auto find_result = std::find(results.begin(),results.end(),index);
                if(find_result != results.end()) {
                    intersection_point_in_cell.push_back(inter_p);
    
                }
            }

            if( v == nullptr){
                return generate_polygon_with_rect_vertexes_and_line(boundary_vertexes_in_cell,intersection_point_in_cell);
            }else{
                Direction direction = e->vertex0() == nullptr ? Direction::in : Direction::out;

                if(!intersection_point_in_cell.empty()){
                    if( direction == Direction::out){
                        polygon_points.emplace_back(v->x(),v->y());

                        polygon_points.emplace_back(intersection_point_in_cell[0].x,intersection_point_in_cell[0].y);
                    }else{
                        polygon_points.emplace_back(intersection_point_in_cell[0].x,intersection_point_in_cell[0].y);

                        polygon_points.emplace_back(v->x(),v->y());
                    }
                }
            }



        }
        e = e->next();
    }while(e != begin);

    auto result = generate_polygon_with_rect_vertexes_and_line(boundary_vertexes_in_cell,polygon_points);
    intersection_polygon(m_boundary,result,polygon_points);

    return polygon_points;
}