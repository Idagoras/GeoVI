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



private:
    geovi::geo::map::GeoMap& gmap;

};

// class GeoMapHandler

void GeoMapHandler::way(const osmium::Way& way){

  //  osmium::WayNodeList& nodes = way.nodes();
  //  osmium::NodeRef* beginNode = const_cast<osmium::NodeRef*>(nodes.begin());
  //  osmium::NodeRef* endNode = const_cast<osmium::NodeRef*>(nodes.end());

    // gmap.addWay();
}

void GeoMapHandler::node(const osmium::Node& node){
   // std::cout << "node id " << node.id() << "latitude = "<< node.location().lat() << " " << "longtitude =" << node.location().lon() << std::endl;
   // gmap.addNode("",node.id());
   TagProcessor tp;
   auto features = tp.getMapFeature(node.tags());
   if( tp.getName(node.tags(),TagProcessor::Language::zh).compare("nameless") != 0){
        std::cout << "node id " << node.id() << " node name is " << tp.getName(node.tags(),TagProcessor::Language::zh) << std::endl;
        for(auto feature : features){
        std::cout << "node's mapfeature is " << get<0>(feature) << " value is " << get<2>(feature) << std::endl;
        }
   }
   GeoMap::GeoNode gnode ;
   gnode.name = tp.getName(node.tags(),TagProcessor::Language::zh);
   gnode.features = tp.getMapFeature(node.tags());
   gnode.id = node.id();
   gnode.loc = GeoMap::Location{node.location().lon(),node.location().lat()};
   gmap.addNode(gnode);
   
  // std::cout << "node id " << node.id() << " node name is " << tp.getName(node.tags(),TagProcessor::Language::zh) << std::endl;
  
  
}

void GeoMapHandler::relation(const osmium::Relation& relation){

}


// class GeoMap

GeoMap::GeoMap(geovi::io::Reader& reader,GeoMapShapeType type,Shape shape):shape_type(type),mshape(shape){
    graph = Graph(0);
    nodes_num = 0;
    ways_num = 0;
    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();

}

bool GeoMap::addNode(GeoNode node){
    nodemap.insert(pair<map_object_id_type,GeoNode>(node.id,node));
    VertexDescriptor v = add_vertex(graph);
    VertexNameMap name_map = get(vertex_name,graph);
    nodes_num ++;
    name_map[v] = std::to_string(nodes_num);
    return true;
}

bool GeoMap::addWay(GeoNode source,GeoNode target){
    ways_num ++;
    return true;
}