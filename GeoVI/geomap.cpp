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
#include "language.h"



using namespace geovi::geo::map;
using namespace geovi::algorithm;
using namespace geovi;
using namespace std;

class TagProcessor {
public:
    using Language = geovi::language::Language;
    static OSMMapFeature getMapFeature(osmium::TagList& tag_list);
    static string getName(osmium::TagList& tag_list,Language language);
    
    
};

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
}

void GeoMapHandler::relation(const osmium::Relation& relation){

}


// class GeoMap
/*
using VertexProperties = property<vertex_name_t,std::string,
                        property<vertex_index_t,int64_t,
                        property<vertex_semantic_sensitivity_t,double,
                        property<vertex_location_t,GeoMap::Location,
                        property<vertex_predecessor_t,vertex_descriptor,
                        property<vertex_distance_t, double>>>>>>;

using EdgeProperties =  property<edge_name_t,std::string,
                        property<edge_index_t,int64_t,
                        property<edge_weight_t,double>>>;

using Graph = adjacency_list<vecS,vecS,directedS,VertexProperties,EdgeProperties>;

using VertexDescriptor = graph_traits<Graph>::vertex_descriptor;
using EdgeDescriptor = graph_traits<Graph>::edge_descriptor;
using VertexIterator = graph_traits<Graph>::vertex_iterator;
using EdgeIterator = graph_traits<Graph>::edge_iterator;

using VertexNameMap = property_map<Graph,vertex_name_t>::type ;
using ConstVertexNameMap = property_map<Graph,vertex_name_t>::const_type ;  
using VertexIndexMap = property_map<Graph,vertex_index_t>::type ;
using ConstVertexIndexMap = property_map<Graph,vertex_index_t>::const_type ;
using VertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::type;
using ConstVertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::const_type;
using VertexLocationMap = property_map<Graph,vertex_location_t>::type;
using ConstVertexLocationMap = property_map<Graph,vertex_location_t>::const_type;


using EdgeNameMap = property_map<Graph,edge_name_t>::type ;
using ConstEdgeNameMap = property_map<Graph,edge_name_t>::const_type ;
using EdgeIndexMap = property_map<Graph,edge_index_t>::type;
using ConstEdgeIndexMap = property_map<Graph,edge_index_t>::const_type;
using EdgeCapacityMap = property_map<Graph,edge_capacity_t>::type;
using ConstEdgeCapacityMap = property_map<Graph,edge_capacity_t>::const_type;
*/
GeoMap::GeoMap(geovi::io::Reader& reader,GeoMapShapeType type,Shape shape):shape_type(type),mshape(shape){
    graph = Graph(0);
    nodes_num = 0;
    ways_num = 0;
    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();
    std::cout<< "node count = " << nodes_num << " " << "way count = " << ways_num << std::endl;
    std::cout << "graph vertex : " << boost::num_vertices(graph) << std::endl;
}

bool GeoMap::addNode(GeoNode node){
    VertexDescriptor v = add_vertex(graph);
    VertexNameMap name_map = get(vertex_name,graph);

   // std::cout<<nodes_num<<std::endl;
    nodes_num ++;
    name_map[v] = std::to_string(nodes_num);
    return true;
}

bool GeoMap::addWay(GeoNode source,GeoNode target){
    ways_num ++;
    return true;
}