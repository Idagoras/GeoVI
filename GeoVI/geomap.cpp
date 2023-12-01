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

Point2 getGridIndex(double max_x,double min_x,
                    double max_y,double min_y,
                    double x,double y,
                    double grid_x,double grid_y){
    return Point2{
            ceil(abs(x-min_x)/grid_x),
            ceil(abs(y-min_y)/grid_y)
    };
}

void init_grid(std::vector<std::vector<GeoMap::GeoGrid>>& grids,
               double max_x,double min_x,
               double max_y,double min_y,
               double grid_x,double grid_y){
    int x_num = static_cast<int>(ceil(abs(max_x-min_x)/grid_x));
    int y_num = static_cast<int>(ceil(abs(max_y-min_y)/grid_y));
    for( int i = 0; i < y_num ; ++ i ){
        std::vector<GeoMap::GeoGrid> grids_x;
        for( int k = 0 ; k < x_num ; ++ k){
            GeoMap::GeoGrid grid;
            grid.index_xy = Point2{static_cast<double>(k),static_cast<double>(i)};
            grids_x.push_back(grid);
        }
        grids.push_back(grids_x);

    }
}

void add_node_index_to_grids(std::vector<std::vector<GeoMap::GeoGrid>>& grids,
                             GeoMap::GeoNode* node_ptr,
                             double max_x,double min_x,
                             double max_y,double min_y,
                             double grid_x,double grid_y){
    auto index_xy = getGridIndex(max_x,min_x,max_y,min_y,node_ptr->utm_xy.x,node_ptr->utm_xy.y,grid_width_1,grid_length_1);
    auto grid_ptr = & grids[static_cast<int>(index_xy.y)][static_cast<int>(index_xy.x)];
    grid_ptr -> nodes.push_back(node_ptr);

}

void add_way_index_to_grids(
        std::vector<std::vector<GeoMap::GeoGrid>>& grids,
        GeoMap::GeoWay* way_ptr,
        double max_x,double min_x,
        double max_y,double min_y,
        double grid_x,double grid_y
        ){
    auto src_node_ptr = &way_ptr ->source;
    auto trg_node_ptr = &way_ptr -> target;
    Point2 src_index_xy = getGridIndex(max_x,min_x,max_y,min_y,src_node_ptr->utm_xy.x,src_node_ptr->utm_xy.y,grid_width_1,grid_length_1);
    Point2 trg_index_xy = getGridIndex(max_x,min_x,max_y,min_y,trg_node_ptr->utm_xy.x,trg_node_ptr->utm_xy.y,grid_width_1,grid_length_1);
    auto src_grid_ptr = & grids[static_cast<int>(src_index_xy.y)][static_cast<int>(src_index_xy.x)];
    auto trg_grid_ptr = & grids[static_cast<int>(trg_index_xy.y)][static_cast<int>(trg_index_xy.x)];
    src_grid_ptr -> ways.push_back(way_ptr);
    trg_grid_ptr -> ways.push_back(way_ptr);
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
            if( name.compare(string(key)) == 0 ){
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
    return string("nameless");
}

class GeoMapHandler : public osmium::handler::Handler{
public:
                

    GeoMapHandler(geovi::geo::map::GeoMap& map):gmap(map){
        osmium::handler::Handler();

    };
    void way(const osmium::Way& way);
    void node(const osmium::Node& node);
    void relation(const osmium::Relation& relation);
    void getGeoNode(const osmium::Node& node,GeoMap::GeoNode& gnode);



private:
    geovi::geo::map::GeoMap& gmap;


};

// class GeoMapHandler

void GeoMapHandler::getGeoNode(const osmium::Node& node,GeoMap::GeoNode& gnode){
    TagProcessor tp;
    auto features = tp.getMapFeature(node.tags());
    if( tp.getName(node.tags(),TagProcessor::Language::zh).compare("nameless") != 0){
       // std::cout << "node id " << node.id() << " node name is " << tp.getName(node.tags(),TagProcessor::Language::zh) << std::endl;
        for(auto feature : features){
       // std::cout << "node's mapfeature is " << get<TagFeild::map_feature>(feature) << " value is " << get<TagFeild::feature_value>(feature) << std::endl;
        }
    }
    gnode.index = gmap.numOfNodes();
    gnode.name = tp.getName(node.tags(),TagProcessor::Language::zh);
    gnode.features = tp.getMapFeature(node.tags());
    gnode.id = node.id();
    gnode.loc = GeoMap::Location{node.location().lon(),node.location().lat()};
    
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
        //std::cout << "node ref in way id " << way.id() << "ref id is " << node.ref() << std::endl;
        //auto geo_node = gmap.getNode(node.ref());
       // std::cout << "node name is " << geo_node -> name << std::endl;
        //std::cout << "latitude is " << geo_node -> loc.latitude << " longitude is " << geo_node -> loc.longitude << std::endl;

    }


    const GeoMap::GeoNode* primary;
    for (int i = 0 ; i < index ; ++ i ){
        if( i == 0 )
            primary = gmap.getNode(refs[i]);
        else {
            auto current = gmap.getNode(refs[i]);
            GeoMap::GeoWay gway = {
                    *primary,
                    *current,
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
    graph = Graph(0);
    nodes_num = 0;
    ways_num = 0;

    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();
    init_grid(m_grids,m_max_x,m_min_x,m_max_y,m_min_y,grid_width_1,grid_length_1);

}

bool GeoMap::addNode(GeoNode node){
    if( !hasNode(node.id)){
        m_node_map.insert(pair<map_object_id_type,GeoNode>(node.id, node));
        m_geo_nodes.push_back(&m_node_map[node.id]);
        addNodeToGraph(node);        
        nodes_num ++;
        node.utm_xy = Point2{ node.loc.latitude,node.loc.longitude};
        s_coordinate_system_converter->convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,node.utm_xy);
        m_max_x = node.utm_xy.x > m_max_x ? node.utm_xy.x : m_max_x;
        m_max_y = node.utm_xy.y > m_max_y ? node.utm_xy.y : m_max_y;
        m_min_x = node.utm_xy.x < m_min_x ? node.utm_xy.x : m_min_x;
        m_min_y = node.utm_xy.y < m_min_y ? node.utm_xy.y : m_min_y;
        return true;
    }
    return false;
}

bool GeoMap::addWay(GeoWay way){
   // std::cout << "sp lat is " << sp.x << " and lon is " << sp.y << std::endl;
   // std::cout << "tp lat is " << tp.x << " and lon is " << tp.y <<std::endl;
    way.capacity = DistanceCalculator::euclidDistance2D(way.source.utm_xy,way.target.utm_xy);
   // std::cout << "way capacity is " << way.capacity << std::endl;
    way.index = ways_num;
    addWayToGraph(way);
    ways_num ++;
    return true;
}

bool GeoMap::hasNode(GeoMap::map_object_id_type node_id){
    return m_node_map.find(node_id) != m_node_map.end();
}

const GeoMap::GeoNode* GeoMap::getNode(GeoMap::map_object_id_type node_id){
    if( hasNode(node_id) ){
        return &(m_node_map.find(node_id)->second);
    }
    return NULL;
}
void GeoMap::addNodeToGraph(GeoNode& node){
    VertexDescriptor v = add_vertex(graph);
    auto name_map = get(vertex_name,graph);
    name_map[v] = node.name;

    auto location_map = get(vertex_location,graph);
    location_map[v] = Point2{node.loc.latitude,node.loc.longitude};
  //  std::cout << "boost graph has " <<graph.m_vertices.size() << " nodes " << std::endl;
}

void GeoMap::addWayToGraph(GeoWay& way){
    int sr_index = way.source.index;
    int tg_index = way.target.index;
    auto sr_vertex_descriptor = vertex(sr_index,graph);
    auto tg_vertex_descriptor = vertex(tg_index,graph);
    EdgeDescriptor e = add_edge(sr_vertex_descriptor,tg_vertex_descriptor,graph).first;
    auto weight_map = get(edge_weight,graph);
    weight_map[e] = way.capacity;
    auto name_map = get(edge_name,graph);
    name_map[e] = way.name;

   //std::cout << "boost graph has " <<graph.m_edges.size() << " ways " << std::endl;
}

std::vector<GeoMap::GeoNode *> GeoMap::getGeoNodes() {
   return m_geo_nodes;
}

std::vector<Point2> GeoMap::getNodes() {
    std::vector<Point2> nodes;
    auto geo_nodes = getGeoNodes();
    for(auto node : geo_nodes){
        nodes.push_back(Point2{node->loc.latitude,node->loc.longitude});
    }
    return nodes;
}