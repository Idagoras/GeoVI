#include "geomap.h"
#include <osmium/handler.hpp>
#include <osmium/visitor.hpp>
#include <iostream>
#include <map>
#include <osmium/osm/way.hpp>
#include <osmium/osm/node_ref_list.hpp>
#include <osmium/osm/node.hpp>
#include <vector>
#include <utility>
#include <tuple>
#include <cmath>
#include "language.h"
#include "convert.h"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <deque>





using namespace geovi::geo::map;
using namespace geovi::algorithm;
using namespace geovi;
using namespace std;


// 1级网格索引
static const double grid_width_1 = 100.0;
static const double grid_length_1 = 100.0;
// 2级网格索引
static const double grid_width_2 = 25.0;
static const double grid_length_2 = 25.0;

static vector<pair<string,OSMMapFeature>> feature_with_name_string = {
    make_pair("amenity",OSMMapFeature::Amenity),
    make_pair("aeroway",OSMMapFeature::Aeroway),
    make_pair("barriers",OSMMapFeature::Barriers),
    make_pair("boundary",OSMMapFeature::Boundary),
    make_pair("building",OSMMapFeature::Building),
    make_pair("craft",OSMMapFeature::Craft),
    make_pair("emergency",OSMMapFeature::Emergency),
    make_pair("geological",OSMMapFeature::Geological),
    make_pair("healthcare",OSMMapFeature::Healthcare),
    make_pair("highway",OSMMapFeature::Highway),
    make_pair("historic",OSMMapFeature::Historic),
    make_pair("landuse",OSMMapFeature::Landuse),
    make_pair("leisure",OSMMapFeature::Leisure),
    make_pair("man_made",OSMMapFeature::Man_made),
    make_pair("military",OSMMapFeature::Military),
    make_pair("natural",OSMMapFeature::Natural),
    make_pair("office",OSMMapFeature::Office),
    make_pair("place",OSMMapFeature::Place),
    make_pair("power",OSMMapFeature::Power),
    make_pair("public_transport",OSMMapFeature::Public_Transport),
    make_pair("railway",OSMMapFeature::Railway),
    make_pair("route",OSMMapFeature::Route),
    make_pair("shop",OSMMapFeature::Shop),
    make_pair("sport",OSMMapFeature::Sport),
    make_pair("telecom",OSMMapFeature::Telecom),
    make_pair("tourism",OSMMapFeature::Tourism),
    make_pair("water",OSMMapFeature::Water),
    make_pair("waterway",OSMMapFeature::Waterway)
};

static std::unique_ptr<CoordinateSystemConverter> s_coordinate_system_converter  =  std::make_unique<CoordinateSystemConverter>(LongitudeBands::band_50);

using GridVertexLocation = enum GridPointLocation {
    left_top = 0,
    left_bottom = 1,
    right_top = 2,
    right_bottom = 3
};

using namespace boost::geometry::strategy::buffer;
using namespace boost::geometry::model;

Point2 getGridIndex(double max_x,double min_x,
                    double max_y,double min_y,
                    double x,double y,
                    double grid_x,double grid_y){
    return Point2{
            floor(abs(y-min_y)/grid_y),
            floor(abs(x-min_x)/grid_x)

    };
}

std::vector<Point2> getUTMCoordinatePoint(double max_x,double min_x,
                                          double max_y,double min_y,
                                          double x,double y,
                                          double grid_x,double grid_y){

    Point2 left_top = Point2{ min_x + y * grid_x,min_y + (x+1) * grid_y };
    Point2 left_bottom = Point2{ min_x + y * grid_x, min_y + x * grid_y };
    Point2 right_top = Point2{ min_x + (y+1) * grid_x, min_y + (x+1) * grid_y};
    Point2 right_bottom = Point2{min_x + (y+1) * grid_x, min_y + x * grid_y };
    return std::vector<Point2>{
            left_top,
            right_top,
            right_bottom,
            left_bottom
    };
}

std::pair<bool,bool> is_within_s_between_c(double grid_index_x,double grid_index_y,double radius,double origin_x,double origin_y,double max_x,double min_x,
                           double max_y,double min_y,double grid_x,double grid_y){
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::polygon<point_xy> polygon;

    const double buffer_distance = radius == 0 ? 0.2 : radius ;
    const int points_per_circle = 36;
    distance_symmetric<double> distance_strategy(buffer_distance);
    join_round join_strategy(points_per_circle);
    end_round end_strategy(points_per_circle);
    point_circle circle_strategy(points_per_circle);
    side_straight side_strategy;

    multi_polygon<polygon> result;
    point_xy origin(origin_x,origin_y);

    boost::geometry::buffer(origin,result,distance_strategy,side_strategy,join_strategy,end_strategy,circle_strategy);

    std::vector<Point2> grid_vertexes = getUTMCoordinatePoint(max_x,min_x,max_y,min_y,grid_index_x,grid_index_y,grid_x,grid_y);
    polygon grid_rect;
    for(auto v : grid_vertexes){
        boost::geometry::append(grid_rect,point_xy(v.x,v.y));
    }

    boost::geometry::append(grid_rect,point_xy(grid_vertexes[0].x,grid_vertexes[0].y));
    return {boost::geometry::within(grid_rect,result.front()) , boost::geometry::within(result.front(),grid_rect)};
}

std::pair<bool,bool> is_within_s_between_p(double grid_index_x,double grid_index_y,double max_x,double min_x,
                           double max_y,double min_y,double grid_x,double grid_y,std::vector<Point2>& polygon_vertexes)
{
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::polygon<point_xy> polygon;

    polygon m_polygon;
    for(auto v : polygon_vertexes){
        boost::geometry::append(m_polygon,point_xy(v.x,v.y));
    }
    boost::geometry::append(m_polygon,point_xy(polygon_vertexes[0].x,polygon_vertexes[0].y));

    std::vector<Point2> grid_vertexes = getUTMCoordinatePoint(max_x,min_x,max_y,min_y,grid_index_x,grid_index_y,grid_x,grid_y);
    polygon grid_rect;
    for(auto v : grid_vertexes){
        boost::geometry::append(grid_rect,point_xy(v.x,v.y));
    }
    boost::geometry::append(grid_rect,point_xy(grid_vertexes[0].x,grid_vertexes[0].y));
    return {boost::geometry::within(grid_rect,m_polygon),boost::geometry::within(m_polygon,grid_rect)};
}

bool is_intersection_s_between_c(double grid_index_x,double grid_index_y,double radius,double origin_x,double origin_y,double max_x,double min_x,
                     double max_y,double min_y,double grid_x,double grid_y){
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::polygon<point_xy> polygon;

    const double buffer_distance = radius == 0 ? 0.2 : radius;
    const int points_per_circle = 36;
    distance_symmetric<double> distance_strategy(buffer_distance);
    join_round join_strategy(points_per_circle);
    end_round end_strategy(points_per_circle);
    point_circle circle_strategy(points_per_circle);
    side_straight side_strategy;

    multi_polygon<polygon> result;
    point_xy origin(origin_x,origin_y);
    boost::geometry::buffer(origin,result,distance_strategy,side_strategy,join_strategy,end_strategy,circle_strategy);

    std::vector<Point2> grid_vertexes = getUTMCoordinatePoint(max_x,min_x,max_y,min_y,grid_index_x,grid_index_y,grid_x,grid_y);
    polygon grid_rect;
    for(auto v : grid_vertexes){
        boost::geometry::append(grid_rect,point_xy(v.x,v.y));
    }
    boost::geometry::append(grid_rect,point_xy(grid_vertexes[0].x,grid_vertexes[0].y));
    std::deque<polygon> intersection_geometries;
    boost::geometry::intersection(grid_rect,result.front(),intersection_geometries);
    if(intersection_geometries.size() == 1){
        return true;
    }
    return false;
}

bool is_interaction_s_between_p(double grid_index_x,double grid_index_y,double max_x,double min_x,
                                double max_y,double min_y,double grid_x,double grid_y,std::vector<Point2>& polygon_vertexes){
    typedef boost::geometry::model::d2::point_xy<double> point_xy;
    typedef boost::geometry::model::polygon<point_xy> polygon;

    polygon m_polygon;
    for(auto v : polygon_vertexes){
        boost::geometry::append(m_polygon,point_xy(v.x,v.y));
    }
    boost::geometry::append(m_polygon,point_xy(polygon_vertexes[0].x,polygon_vertexes[0].y));

    std::vector<Point2> grid_vertexes = getUTMCoordinatePoint(max_x,min_x,max_y,min_y,grid_index_x,grid_index_y,grid_x,grid_y);
    polygon grid_rect;
    for(auto v : grid_vertexes){
        boost::geometry::append(grid_rect,point_xy(v.x,v.y));
    }
    boost::geometry::append(grid_rect,point_xy(grid_vertexes[0].x,grid_vertexes[0].y));
    std::deque<polygon> intersection_geometries;
    boost::geometry::intersection(grid_rect,m_polygon,intersection_geometries);
    if(intersection_geometries.size() == 1){
        return true;
    }
    return false;
}


std::pair<LongitudeBands,LongitudeBands> get_UTM_longitude_bands_range(double max_lon,double min_lon,double min_lat,double max_lat){
    return {s_coordinate_system_converter->utm_identity(min_lon,max_lat),
            s_coordinate_system_converter->utm_identity(max_lon,max_lat)};
}




void init_grid(std::vector<std::vector<GeoMap::GeoGrid>>& grids,
                double max_x,double min_x,
                double max_y,double min_y,
                double grid_x,double grid_y){
    int x_num = static_cast<int>(ceil(abs(max_x-min_x)/grid_x));
    int y_num = static_cast<int>(ceil(abs(max_y-min_y)/grid_y));
    std::cout << "grids num = " << x_num * y_num << std::endl;
    for( int i = 0; i < y_num ; ++ i ){
        std::vector<GeoMap::GeoGrid> grids_x;
        for( int k = 0 ; k < x_num ; ++ k){
            GeoMap::GeoGrid grid;
            grid.index_xy = Point2{static_cast<double>(i),static_cast<double>(k)};
            grids_x.push_back(grid);
        }
        grids.push_back(grids_x);

    }
    /*
    for(auto v_x : grids){
        for(auto x : v_x){
            std::cout << "grid.index x = " << x.index_xy.x << " and index y = " << x.index_xy.y << std::endl;
        }
    }
     */
}

void add_node_index_to_grids(std::vector<std::vector<GeoMap::GeoGrid>>& grids,
                            GeoMap::GeoNode* node_ptr,
                            double max_x,double min_x,
                            double max_y,double min_y,
                            double grid_x,double grid_y){
    auto index_xy = getGridIndex(max_x,min_x,max_y,min_y,node_ptr->utm_xy.x,node_ptr->utm_xy.y,grid_x,grid_y);
    auto grid_ptr = & grids[static_cast<int>(index_xy.x)][static_cast<int>(index_xy.y)];
    grid_ptr -> nodes.push_back(node_ptr);

}

void add_way_index_to_grids(
        std::vector<std::vector<GeoMap::GeoGrid>>& grids,
        GeoMap::GeoWay* way_ptr,
        double max_x,double min_x,
        double max_y,double min_y,
        double grid_x,double grid_y
        ){
    auto src_node_ptr = way_ptr -> source;
    auto trg_node_ptr = way_ptr -> target;
    Point2 src_index_xy = getGridIndex(max_x,min_x,max_y,min_y,src_node_ptr->utm_xy.x,src_node_ptr->utm_xy.y,grid_x,grid_y);
    Point2 trg_index_xy = getGridIndex(max_x,min_x,max_y,min_y,trg_node_ptr->utm_xy.x,trg_node_ptr->utm_xy.y,grid_x,grid_y);
    auto src_grid_ptr = & grids[static_cast<int>(src_index_xy.x)][static_cast<int>(src_index_xy.y)];
    auto trg_grid_ptr = & grids[static_cast<int>(trg_index_xy.x)][static_cast<int>(trg_index_xy.y)];
    src_grid_ptr -> ways.push_back(way_ptr);
    trg_grid_ptr -> ways.push_back(way_ptr);
}

void fill_grids(std::vector<std::vector<GeoMap::GeoGrid>>& grids,
                std::vector<GeoMap::GeoWay*>& way_ptrs,
                std::vector<GeoMap::GeoNode*>& node_ptrs,
                double max_x,double min_x,
                double max_y,double min_y,
                double grid_x,double grid_y){

    for(auto way_ptr : way_ptrs){
        add_way_index_to_grids(grids,way_ptr,max_x,min_x,max_y,min_y,grid_x,grid_y);
    }
    std::cout << "have added way to grids " << std::endl;
    for(auto node_ptr : node_ptrs){
        add_node_index_to_grids(grids,node_ptr,max_x,min_x,max_y,min_y,grid_x,grid_y);
    }
    std::cout << "have added node to grids " << std::endl;
}

using TagField = enum TagField {
    name = 0,
    map_feature = 1,
    feature_value = 2
};

class TagProcessor {
public:
    using Language = geovi::language::Language;
    using FeatureTupleArray = vector<std::tuple<string,OSMMapFeature,string>>;
    FeatureTupleArray getMapFeature(const osmium::TagList& tag_list);
    string getName(const osmium::TagList& tag_list,Language language);
    
    
};

// class TagProcessor

TagProcessor::FeatureTupleArray TagProcessor::getMapFeature(const osmium::TagList& tag_list){
    FeatureTupleArray arr;
    for(auto it = tag_list.begin(); it != tag_list.end(); it ++ ){
        auto key = it -> key();
        bool exist = false;
        for(auto it_f = feature_with_name_string.begin() ; it_f != feature_with_name_string.end(); it_f ++){
            string name = it_f->first;
            if( name == string(key) ){
                string value = it -> value();
                arr.push_back(make_tuple(it_f->first,it_f->second,value));
                break;
            }
        }
    }
    return arr;

}

string TagProcessor::getName(const osmium::TagList& tag_list,Language language){
    if( tag_list.has_key("name") ){
        return tag_list["name"];
    }
    return {"nameless"};
}

class GeoMapHandler : public osmium::handler::Handler{
public:
                

    explicit GeoMapHandler(geovi::geo::map::GeoMap& map):gmap(map){
    };
    void way(const osmium::Way& way);
    void node(const osmium::Node& node);
    void relation(const osmium::Relation& relation);
    void getGeoNode(const osmium::Node& node,GeoMap::GeoNode& g_node);



private:
    geovi::geo::map::GeoMap& gmap;



};

// class GeoMapHandler

void GeoMapHandler::getGeoNode(const osmium::Node& node,GeoMap::GeoNode& g_node){
    TagProcessor tp;
    auto features = tp.getMapFeature(node.tags());
    g_node.index = gmap.numOfNodes();
    g_node.name = tp.getName(node.tags(),TagProcessor::Language::zh);
    g_node.features = tp.getMapFeature(node.tags());
    g_node.id = node.id();
    g_node.loc = GeoMap::Location{node.location().lon(),node.location().lat()};
    
}

void GeoMapHandler::way(const osmium::Way& way){

    //std::cout << "way's id is " << way.id() << std::endl;
    TagProcessor tp;
    auto features = tp.getMapFeature(way.tags());
    auto name = tp.getName(way.tags(),TagProcessor::Language::zh);
    const osmium::WayNodeList& nodes = way.nodes();
    osmium::object_id_type refs[nodes.size()];
    int index = 0;                                                                                                                                                        
    for(auto it = nodes.begin(); it != nodes.end(); ++ it ){
        auto node = *it;
        if( !gmap.hasNode(node.ref())){
            std::cout << "node is not exist in nodes set" <<  std::endl;
            std::cout << "fail to add way" << std::endl;
            return;
        }
        refs[index] = node.ref();
        ++ index;

    }


    const GeoMap::GeoNode* primary;
    for (int i = 0 ; i < index ; ++ i ){
        if( i == 0 )
            primary = gmap.getNode(refs[i]);
        else {
            auto current = gmap.getNode(refs[i]);
            GeoMap::GeoWay gway = {
                    primary,
                    current,
                    way.id(),
                    0,
                    -1,
                    name,
                    features
            };
            gmap.addWay(gway);
            primary = current;
        }
    }
}

void GeoMapHandler::node(const osmium::Node& node){
   // std::cout << "node id " << node.id() << "latitude = "<< node.location().lat() << " " << "longitude =" << node.location().lon() << std::endl;
   // gmap.addNode("",node.id());

    
    if ( !gmap.hasNode(node.id())){
        GeoMap::GeoNode gnode;
        getGeoNode(node,gnode);
        gmap.addNode(gnode);
    
    }

  // std::cout << "node id " << node.id() << " node name is " << tp.getName(node.tags(),TagProcessor::Language::zh) << std::endl;


}

void GeoMapHandler::relation(const osmium::Relation& relation){

}



// class GeoMap

GeoMap::GeoMap(geovi::io::OSMReader& reader,GeoMapShapeType type,Shape shape):shape_type(type),mshape(shape){
    m_graph = Graph(0);
    m_nodes_num = 0;
    m_ways_num = 0;
    m_relations_num = 0;

    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();
    init_grid(m_grids,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);
    fill_grids(m_grids,m_geo_ways,m_geo_nodes,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);

}

GeoMap::GeoMap(geovi::io::OSMReader& reader,const std::shared_ptr<MapFeatureFilter>& filter):m_filter(filter){
    m_graph = Graph(0);
    m_nodes_num = 0;
    m_ways_num = 0;
    m_relations_num = 0;
    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();
    auto range_pair = get_UTM_longitude_bands_range(m_max_lon,m_min_lon,m_min_lat,m_max_lat);
    std::cout << "min bands num = " << range_pair.first << " and max bands num =" << range_pair.second << std::endl << std::endl;
    init_grid(m_grids,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);
    fill_grids(m_grids,m_geo_ways,m_geo_nodes,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);
    std::cout << "have filled grids in map" <<std::endl;
}

bool GeoMap::addNode(GeoNode node){
    if( !hasNode(node.id)){

        m_node_map.insert(pair<map_object_id_type,GeoNode>(node.id, node));
        m_geo_nodes.push_back(&m_node_map[node.id]);
        m_add_node_to_graph(node);
        m_nodes_num ++;
        m_node_map[node.id].utm_xy = Point2{ node.loc.latitude,node.loc.longitude};
        s_coordinate_system_converter->convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,m_node_map[node.id].utm_xy);
        m_max_x = m_node_map[node.id].utm_xy.x > m_max_x ? m_node_map[node.id].utm_xy.x : m_max_x;
        m_max_y = m_node_map[node.id].utm_xy.y > m_max_y ? m_node_map[node.id].utm_xy.y : m_max_y;
        m_min_x = m_node_map[node.id].utm_xy.x < m_min_x ? m_node_map[node.id].utm_xy.x : m_min_x;
        m_min_y = m_node_map[node.id].utm_xy.y < m_min_y ? m_node_map[node.id].utm_xy.y : m_min_y;
        m_max_lat = m_node_map[node.id].loc.latitude > m_max_lat ? m_node_map[node.id].loc.latitude : m_max_lat ;
        m_max_lon = m_node_map[node.id].loc.longitude > m_max_lon ? m_node_map[node.id].loc.longitude : m_max_lon ;
        m_min_lat = m_node_map[node.id].loc.latitude < m_min_lat ? m_node_map[node.id].loc.latitude : m_min_lat ;
        m_min_lon = m_node_map[node.id].loc.longitude < m_min_lon ? m_node_map[node.id].loc.longitude : m_min_lon ;

        if( !m_filter.expired() ){
            auto filter = m_filter.lock();
            for(auto& map_feature_tuple : m_node_map[node.id].features ){
                filter ->node_filter(get<TagField::map_feature>(map_feature_tuple),get<TagField::feature_value>(map_feature_tuple).c_str(),&m_node_map[node.id]);
            }
        }
        return true;
    }
    return false;
}

bool GeoMap::addWay(GeoWay way){
    if(!hasWay(way.id)) {
       // std::cout << "add way id = " << way.id << std::endl;
        m_way_map.insert({way.id,way});
        m_geo_ways.push_back(&m_way_map[way.id]);
        way.capacity = DistanceCalculator::euclidDistance2D(way.source->utm_xy, way.target->utm_xy);
        way.index = m_ways_num;
        m_add_way_to_graph(way);
        m_ways_num++;

        if( !m_filter.expired() ){
            auto filter = m_filter.lock();
            for(auto& map_feature_tuple : way.features ){
                filter ->way_filter(get<TagField::map_feature>(map_feature_tuple),get<TagField::feature_value>(map_feature_tuple).c_str(),&m_way_map[way.id]);
            }
        }
        return true;
    }
    return false;
}

bool GeoMap::hasNode(GeoMap::map_object_id_type node_id){
    return m_node_map.find(node_id) != m_node_map.end();
}

bool GeoMap::hasWay(geovi::geo::map::GeoMap::map_object_id_type way_id) {
    return m_way_map.find(way_id) != m_way_map.end();
}

const GeoMap::GeoNode* GeoMap::getNode(GeoMap::map_object_id_type node_id){
    if( hasNode(node_id) ){
        return &(m_node_map.find(node_id)->second);
    }
    return nullptr;
}

void GeoMap::m_add_node_to_graph(GeoNode& node){
    VertexDescriptor v = add_vertex(m_graph);
    auto name_map = get(vertex_name, m_graph);
    name_map[v] = node.name;

    auto location_map = get(vertex_location, m_graph);
    location_map[v] = Point2{node.loc.latitude,node.loc.longitude};

}

void GeoMap::m_add_way_to_graph(GeoWay& way){
    auto sr_index = way.source->index;
    auto tg_index = way.target->index;
    auto sr_vertex_descriptor = vertex(sr_index, m_graph);
    auto tg_vertex_descriptor = vertex(tg_index, m_graph);
    EdgeDescriptor e = add_edge(sr_vertex_descriptor, tg_vertex_descriptor, m_graph).first;
    auto weight_map = get(edge_weight, m_graph);
    weight_map[e] = way.capacity;
    auto name_map = get(edge_name, m_graph);
    name_map[e] = way.name;

}

std::vector<GeoMap::GeoNode *> GeoMap::getGeoNodes() {
    return m_geo_nodes;
}

std::vector<Point2> GeoMap::getUTMNodesCoordinate() {
    std::vector<Point2> nodes;
    auto geo_nodes = getGeoNodes();
    for(auto node : geo_nodes){
        nodes.emplace_back(node->loc.latitude,node->loc.longitude);
    }
    return nodes;
}

std::vector<const GeoMap::GeoNode *> GeoMap::find(int radius_metre, double utm_x, double utm_y) {
    Point2 origin = {utm_x,utm_y};
    std::vector<const GeoMap::GeoNode *> find_node_ptrs;
    std::vector<Point2> fill_grids,un_fill_grids;
    auto right_grid_index = getGridIndex(m_max_x,m_min_x,m_max_y,m_min_y,utm_x + radius_metre,utm_y,grid_width_1,grid_length_1);
    right_grid_index = utm_x + radius_metre > m_max_x ? getGridIndex(m_max_x,m_min_x,m_max_y,m_min_y,m_max_x,utm_y,grid_width_1,grid_length_1) : right_grid_index ;
    auto left_grid_index = getGridIndex(m_max_x,m_min_x,m_max_y,m_min_y,utm_x - radius_metre,utm_y,grid_width_1,grid_length_1);
    left_grid_index = utm_x - radius_metre < m_min_x ? getGridIndex(m_max_x,m_min_x,m_max_y,m_min_y,m_min_x,utm_y,grid_width_1,grid_length_1) : left_grid_index ;
    auto m_bottom_grid_index = getGridIndex(m_max_x, m_min_x, m_max_y, m_min_y, utm_x , utm_y - radius_metre, grid_width_1, grid_length_1);
    m_bottom_grid_index = utm_y - radius_metre < m_min_y ? getGridIndex(m_max_x, m_min_x, m_max_y, m_min_y, utm_x, m_min_y, grid_width_1, grid_length_1) : m_bottom_grid_index ;
    auto m_top_grid_index = getGridIndex(m_max_x, m_min_x, m_max_y, m_min_y, utm_x, utm_y + radius_metre, grid_width_1, grid_length_1);
    m_top_grid_index = utm_y + radius_metre > m_max_y ? getGridIndex(m_max_x, m_min_x, m_max_y, m_min_y, utm_x, m_max_y, grid_width_1, grid_length_1) : m_top_grid_index ;


    for( int i = static_cast<int>(left_grid_index.y); i <= static_cast<int>(right_grid_index.y) ; ++ i ){
        for (int j = static_cast<int>(m_bottom_grid_index.x) ; j <= static_cast<int>(m_top_grid_index.x) ; ++ j ){
            auto is_within_pair = is_within_s_between_c(j,i,radius_metre,utm_x,utm_y,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);
            if( is_within_pair.second){
                un_fill_grids.emplace_back(static_cast<double>(j), static_cast<double>(i));
                break;
            }else{
                if( is_within_pair.first ){
                    fill_grids.emplace_back(static_cast<double>(j), static_cast<double>(i));
                }else{
                    bool is_intersection = is_intersection_s_between_c(j,i,radius_metre,utm_x,utm_y,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);
                    if ( is_intersection )
                        un_fill_grids.emplace_back(static_cast<double>(j), static_cast<double>(i));
                }
            }

        }
    }

    for(auto grid : fill_grids){
        auto grid_ptr = & m_grids[static_cast<int>(grid.x)][static_cast<int>(grid.y)];
        find_node_ptrs.insert(find_node_ptrs.end(),grid_ptr -> nodes.begin(),grid_ptr -> nodes.end());
    }
    for(auto grid : un_fill_grids){
        auto grid_ptr = & m_grids[static_cast<int>(grid.x)][static_cast<int>(grid.y)];
        for( auto node : grid_ptr -> nodes){
            if ( DistanceCalculator::euclidDistance2D(node->utm_xy,origin) <= radius_metre )
                find_node_ptrs.push_back(node);
        }
    }
    return find_node_ptrs;

}

std::vector<double> GeoMap::shortestPathsDistance(double utm_x, double utm_y) {

    std::vector<double> results;

    auto grid_index_xy = getGridIndex(m_max_x,m_min_x,m_max_y,m_min_y,utm_x,utm_y,grid_width_1,grid_length_1);
    auto grid_ptr = &m_grids[static_cast<int>(grid_index_xy.x)][static_cast<int>(grid_index_xy.y)];
    double min_distance = std::numeric_limits<double>::max();
    int64_t index = -1;
    for(auto node : grid_ptr->nodes){
        double distance = DistanceCalculator::euclidDistance2D(node->utm_xy,{utm_x,utm_y});
        if(distance < min_distance){
            min_distance = distance ;
            index = node -> index;
        }
    }
    VertexDescriptor origin_vertex = vertex(index, m_graph);
    VertexPredecessorMap p = get(vertex_predecessor, m_graph);
    VertexDistanceMap d = get(vertex_distance, m_graph);
    dijkstra_shortest_paths(m_graph, origin_vertex, predecessor_map(p).distance_map(d));
    for(auto node_ptr : m_geo_nodes){
        auto v = vertex(node_ptr->index, m_graph);
        std::cout << "origin index " << index << " node index "<< node_ptr -> index <<" distance = " << d[v] << std::endl;
        results.push_back(d[v]+min_distance);
    }

    return results;
}



using Highway_value = enum highway_value {
    bus_stop = 0,
    crossing ,
    cyclist_waiting_aid ,
    elevator ,
    emergency_bay,
    emergency_access_point,
    give_way,
    phone,
    ladder,
    milestone,
    mini_roundabout,
    motorway_junction,
    passing_place,
    platform,
    rest_area,
    services,
    speed_camera,
    stop,
    street_lamp,
    toll_gantry,
    traffic_signals,
    traffic_mirror,
    trailhead,
    turning_circle,
    turning_loop

};

#define FEATURE_VALUE_PAIR(VALUE) {#VALUE,Highway_value::VALUE}

static std::map<std::string,Highway_value> highway_feature_value_string_to_highway_value_enum_map{
    FEATURE_VALUE_PAIR(bus_stop),
    FEATURE_VALUE_PAIR(crossing),
    FEATURE_VALUE_PAIR(cyclist_waiting_aid),
    FEATURE_VALUE_PAIR(elevator),
    FEATURE_VALUE_PAIR(emergency_bay),
    FEATURE_VALUE_PAIR(emergency_access_point),
    FEATURE_VALUE_PAIR(give_way),
    FEATURE_VALUE_PAIR(phone),
    FEATURE_VALUE_PAIR(ladder),
    FEATURE_VALUE_PAIR(milestone),
    FEATURE_VALUE_PAIR(mini_roundabout),
    FEATURE_VALUE_PAIR(motorway_junction),
    FEATURE_VALUE_PAIR(passing_place),
    FEATURE_VALUE_PAIR(platform),
    FEATURE_VALUE_PAIR(rest_area),
    FEATURE_VALUE_PAIR(services),
    FEATURE_VALUE_PAIR(speed_camera),
    FEATURE_VALUE_PAIR(stop),
    FEATURE_VALUE_PAIR(street_lamp),
    FEATURE_VALUE_PAIR(toll_gantry),
    FEATURE_VALUE_PAIR(traffic_signals),
    FEATURE_VALUE_PAIR(traffic_mirror),
    FEATURE_VALUE_PAIR(trailhead),
    FEATURE_VALUE_PAIR(turning_circle),
    FEATURE_VALUE_PAIR(turning_loop)
};

void MapFeatureFilter::node_filter(geovi::geo::map::OSMMapFeature feature, const char *feature_value,
                                   const GeoMap::GeoNode *node) {

}

void MapFeatureFilter::way_filter(geovi::geo::map::OSMMapFeature feature, const char *feature_value,
                                  const GeoMap::GeoWay *node) {

}

void CrossingFilter::node_filter(OSMMapFeature feature,const char * feature_value,const GeoMap::GeoNode* node){
    std::vector<int64_t> feature_value_nums(25);
    if( feature == OSMMapFeature::Highway ){
        if(highway_feature_value_string_to_highway_value_enum_map.find(std::string(feature_value)) != highway_feature_value_string_to_highway_value_enum_map.end()){
            int64_t index = highway_feature_value_string_to_highway_value_enum_map[std::string(feature_value)];
            ++ feature_value_nums[index] ;
            if(std::strcmp(feature_value,"crossing") == 0){
                m_crossing_nodes.push_back(node);
            }
        }

    }
}

std::vector<Point2> CrossingFilter::crossing_utm_xy_points(){
    std::vector<Point2> utm_xy_points ;
    for(auto node_ptr : m_crossing_nodes){
        utm_xy_points.emplace_back(node_ptr->utm_xy.x,node_ptr->utm_xy.y);
    }
    return utm_xy_points ;
}

