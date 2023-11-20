#include "geomap.h"
#include <osmium/handler.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/visitor.hpp>
#include <iostream>
#include <map>
#include <osmium/osm/way.hpp>
#include <osmium/osm/node_ref_list.hpp>
#include <osmium/osm/node.hpp>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include "language.h"



using namespace geovi::geo::map;
using namespace geovi::algorithm;
using namespace geovi;
using namespace std;


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
    make_pair("millitary",OSMMapFeature::Millitary),
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


using TagField = enum TagFeild {
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
        std::cout << "node id " << node.id() << " node name is " << tp.getName(node.tags(),TagProcessor::Language::zh) << std::endl;
        for(auto feature : features){
        std::cout << "node's mapfeature is " << get<TagFeild::map_feature>(feature) << " value is " << get<TagFeild::feature_value>(feature) << std::endl;
        }
    }
    gnode.index = gmap.numOfNodes();
    gnode.name = tp.getName(node.tags(),TagProcessor::Language::zh);
    gnode.features = tp.getMapFeature(node.tags());
    gnode.id = node.id();
    gnode.loc = GeoMap::Location{node.location().lon(),node.location().lat()};
    
}

void GeoMapHandler::way(const osmium::Way& way){
    std::cout << "way's id is " << way.id() << std::endl;
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
    }
    const GeoMap::GeoNode* primary;
    for (int i = 0 ; i < index ; ++ i ){
        if(primary == NULL)
            primary = gmap.getNode(refs[i]);
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

void GeoMapHandler::node(const osmium::Node& node){
   // std::cout << "node id " << node.id() << "latitude = "<< node.location().lat() << " " << "longtitude =" << node.location().lon() << std::endl;
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

}

bool GeoMap::addNode(GeoNode node){
    if( !hasNode(node.id)){
        nodemap.insert(pair<map_object_id_type,GeoNode>(node.id,node));
        addNodeToGraph(node);        
        nodes_num ++;                                                                                          
        return true;
    }
    return false;
}

bool GeoMap::addWay(GeoWay way){
    Point2 sp = {way.source.loc.latitude,way.source.loc.longitude};
    Point2 tp = {way.target.loc.latitude,way.target.loc.longitude};
    CoordinateSystemConverter converter(LongitudeBands::band_31);
    converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,sp);
    converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,tp);
    way.capacity = DistanceCaculator::euclidDistance2D(sp,tp);
    way.index = ways_num;
    addWayToGraph(way);
    ways_num ++;
    return true;
}

bool GeoMap::hasNode(GeoMap::map_object_id_type node_id){
    return nodemap.find(node_id) != nodemap.end();
}

const GeoMap::GeoNode* GeoMap::getNode(GeoMap::map_object_id_type node_id){
    if( hasNode(node_id) ){
        return &(nodemap.find(node_id)->second);
    }
    return NULL;
}
void GeoMap::addNodeToGraph(GeoNode& node){
    VertexDescriptor v = add_vertex(graph);
    VertexNameMap name_map = get(vertex_name,graph);
    name_map[v] = node.name;

    VertexLocationMap location_map = get(vertex_location,graph);
    location_map[v] = Point2{node.loc.latitude,node.loc.longitude};

}

void GeoMap::addWayToGraph(GeoWay& way){
    int sr_index = way.source.index;
    int tg_index = way.target.index;
    auto sr_vertex_descriptor = vertex(sr_index,graph);
    auto tg_vertex_descriptor = vertex(tg_index,graph);
    EdgeDescriptor e = add_edge(sr_vertex_descriptor,tg_vertex_descriptor,graph).first;
    EdgeWeightMap weight_map = get(edge_weight,graph);
    weight_map[e] = way.capacity;
    EdgeNameMap name_map = get(edge_name,graph);
    name_map[e] = way.name;
}