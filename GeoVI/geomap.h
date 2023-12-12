#ifndef GEOVI_GEOMAP_H
#define GEOVI_GEOMAP_H

#include "reader.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/polygon/voronoi.hpp>
#include <vector>
#include <string>
#include "convert.h"
#include <map>
#include <memory>
#include <limits>


namespace boost{
    enum vertex_semantic_sensitivity_t{
        vertex_semantic_sensitivity = 2042
    };
    BOOST_INSTALL_PROPERTY(vertex,semantic_sensitivity);

    enum vertex_location_t{
        vertex_location = 1111
    };
    BOOST_INSTALL_PROPERTY(vertex,location);

}
namespace geovi
{
    namespace algorithm{
        using namespace boost;

        using VertexDescriptor = adjacency_list_traits<vecS,vecS,undirectedS>::vertex_descriptor;
        using EdgeDescriptor = adjacency_list_traits<vecS,vecS,undirectedS>::edge_descriptor;


        using VertexProperties = property<vertex_name_t,std::string,
                        property<vertex_index_t,int64_t,
                        property<vertex_semantic_sensitivity_t,double,
                        property<vertex_location_t,Point2,
                        property<vertex_predecessor_t, VertexDescriptor,
                        property<vertex_distance_t, double>>>>>>;

        using EdgeProperties =  property<edge_name_t,std::string,
                                property<edge_index_t,int64_t,
                                property<edge_weight_t,double>>>;

        using Graph = adjacency_list<vecS,vecS,undirectedS,VertexProperties,EdgeProperties>;



        using VertexNameMap = property_map<Graph,vertex_name_t>::type ;
        using ConstVertexNameMap = property_map<Graph,vertex_name_t>::const_type ;  
        using VertexIndexMap = property_map<Graph,vertex_index_t>::type ;
        using ConstVertexIndexMap = property_map<Graph,vertex_index_t>::const_type ;
        using VertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::type;
        using ConstVertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::const_type;
        using VertexLocationMap = property_map<Graph,vertex_location_t>::type;
        using ConstVertexLocationMap = property_map<Graph,vertex_location_t>::const_type;
        using VertexDistanceMap = property_map<Graph,vertex_distance_t>::type;
        using ConstVertexDistanceMap = property_map<Graph,vertex_distance_t>::const_type;
        using VertexPredecessorMap = property_map<Graph,vertex_predecessor_t>::type;
        using ConstVertexPredecessorMap = property_map<Graph,vertex_predecessor_t>::const_type;


        using EdgeNameMap = property_map<Graph,edge_name_t>::type ;
        using ConstEdgeNameMap = property_map<Graph,edge_name_t>::const_type ;
        using EdgeIndexMap = property_map<Graph,edge_index_t>::type;
        using ConstEdgeIndexMap = property_map<Graph,edge_index_t>::const_type;
        using EdgeWeightMap = property_map<Graph,edge_weight_t>::type;
        using ConstEdgeWeightMap = property_map<Graph,edge_weight_t>::const_type;
    }
    namespace geo
    {
        namespace map{

                

                typedef enum{
                    Amenity, // Used to map facilities used by visitors and residents.
                    Aeroway, // This is used to tag different forms of transportation for people or goods by using aerial wires. 
                    Barriers, // These are used to describe barriers and obstacles that are usually involved by traveling.
                    Boundary, // These are used to describe administrative and other boundaries.
                    Building, // This is used to identify individual buildings or groups of connected buildings.
                    Craft, // This is used as a place that produces or processes customised goods.
                    Emergency, // This is used to describe the location of emergency facilities and equipment. 
                    Geological, // This is used to describe the geological makeup of an area.
                    Healthcare, // Healthcare features.
                    Highway, // This is used to describe roads and footpaths.
                    Historic, // This is used to describe various historic places. 
                    Landuse, // This is used to describe the purpose for which an area of land is being used. 
                    Leisure, // This is used to tag leisure and sports facilities.
                    Man_made, // A tag for identifying man made (artificial) structures that are added to the landscape.
                    Military, // This is used for facilities and on land used by the military.
                    Natural, // This is used to describe natural and physical land features.
                    Office, // An office is a place of business where administrative or professional work is carried out.
                    Place, // This is used mainly to give details about settlements. 
                    Power, // These are used to map electrical power generation and distributions systems. 
                    Public_Transport, // This is used for features related to public transport.
                    Railway, // This tag includes all kinds of railways ranging from heavily used mainline railways to an abandoned rail line.
                    Route, // This is used to describe routes of all different kinds.
                    Shop, // The shop tag is used as a place of business that has stocked goods for sale or sells services. 
                    Sport, // This is used to provide information about which sports are placed on a facility such as a pitch or stadium.
                    Telecom, // These are used to map telecommunication systems.
                    Tourism, // This is used to map places and things of specific interest to tourists.
                    Water, // This is used to describe type of water body and is only used together with natural=water.
                    Waterway, // This is used to described different types of waterways. When mapping the way of a river, stream, drain, canal, etc. these need to be aligned in the direction of the water flow. 
                    None,
                    
                } OSMMapFeature; // Standard OSM Map Features

            class MapFeatureFilter;

            class GeoMap{
            public:
                
                typedef int64_t map_object_id_type ;
            
                typedef enum {
                    circle,
                    rectangle,
                } GeoMapShapeType;
                
                typedef struct{
                    // longitude
                        double longitude;
                    // latitude
                        double latitude;

                } Location;

                typedef struct 
                {
                    // 中心点
                    Location center;
                } Shape;
                

                typedef struct Circle : public Shape
                {
                    // 半径
                    double radius;
                    // 中心点
                    Location center;

                } Circle;

                typedef struct Rectangle : public Shape{
                    // 长
                    double length;
                    // 宽
                    double width;
                    // 中心点
                    Location center;
                    // 左上角的顶点
                    Location top_left_vertex;
                } Rectangle;
                
                typedef struct GeoNode{
                    std::string name;
                    map_object_id_type id;
                    int64_t index = -1;
                    double semantic_sensitivity = -1;
                    Location loc;
                    Point2 utm_xy;
                    std::vector<std::tuple<std::string,OSMMapFeature,std::string>> features;
                } GeoNode;

                typedef struct GeoWay{
                    const GeoNode* source;
                    const GeoNode* target;
                    map_object_id_type id;
                    double capacity = 0;
                    int64_t index = -1;
                    std::string name;
                    std::vector<std::tuple<std::string,OSMMapFeature,std::string>> features;
                } GeoWay;

                typedef struct GeoGrid{
                    std::vector<const GeoNode*> nodes;
                    std::vector<const GeoWay*> ways;
                    Point2 index_xy ;
                } GeoGrid;

                typedef std::map<map_object_id_type,GeoNode> NodeMap;
                typedef std::map<map_object_id_type,GeoWay> WayMap;
                GeoMapShapeType shape_type;
                Shape mshape;
                GeoMap(geovi::io::OSMReader& reader,GeoMapShapeType type,Shape shape);
                GeoMap(geovi::io::OSMReader& reader,const std::shared_ptr<MapFeatureFilter>& filter);

                inline int64_t numOfNodes() const{

                    return m_nodes_num;
                }

                inline int64_t numOfWays() const{
                    return m_ways_num;
                }

                std::vector<const GeoNode* > find(int radius_metre,double utm_x,double utm_y);

                bool addNode(GeoNode node);
                bool addWay(GeoWay way);
                bool hasNode(map_object_id_type node_id);
                bool hasWay(map_object_id_type way_id);
                const GeoNode* getNode(map_object_id_type node_id);
                // 返回地图中所有点的utm坐标集合
                std::vector<Point2> getUTMNodesCoordinate();
                // 返回地图中所有顶点的信息集合
                std::vector<GeoNode*> getGeoNodes();
                std::vector<double> shortestPathsDistance(double utm_x,double utm_y);

                std::vector<Point2> utm_boundary();



               // GeoMap POISOfAnAreaWithCenter(Location loc);
               // GeoMap POISOfAnAreaWithCenterObjectID(map_object_id_type id);


            private:
                geovi::algorithm::Graph m_graph;
                NodeMap m_node_map;
                std::vector<GeoNode*> m_geo_nodes;
                WayMap m_way_map;
                std::vector<GeoWay*> m_geo_ways;
                int m_nodes_num;
                int m_ways_num;
                int m_relations_num;
                double m_max_x = std::numeric_limits<double>::min();
                double m_max_y = std::numeric_limits<double>::min();
                double m_min_x = std::numeric_limits<double>::max();
                double m_min_y = std::numeric_limits<double>::max();
                double m_max_lat = std::numeric_limits<double>::min();
                double m_max_lon = std::numeric_limits<double>::min();
                double m_min_lat = std::numeric_limits<double>::max();
                double m_min_lon = std::numeric_limits<double>::max();
                std::weak_ptr<MapFeatureFilter> m_filter;
                std::vector<std::vector<GeoGrid>> m_grids;
                void m_add_node_to_graph(GeoNode& node);
                void m_add_way_to_graph(GeoWay& way);
            };   

        class  MapFeatureFilter {
        public:
            virtual void node_filter(OSMMapFeature feature,const char * feature_value,const GeoMap::GeoNode* node);
            virtual void way_filter(OSMMapFeature feature,const char * feature_value,const GeoMap::GeoWay* node);
        };

        class CrossingFilter : public MapFeatureFilter {
        public:
            void node_filter(OSMMapFeature feature,const char * feature_value,const GeoMap::GeoNode* node) override;
            std::vector<Point2> crossing_utm_xy_points();
            inline std::vector<const GeoMap::GeoNode*> crossing_geo_nodes(){
                return m_crossing_nodes;
            }
        private:
            std::vector<const GeoMap::GeoNode*> m_crossing_nodes; 

        };


        class GeoMapVoronoiDiagramAdaptor {
        public:
            using voronoi_diagram = boost::polygon::voronoi_diagram<double>;
            using voronoi_cell_container_type = boost::polygon::voronoi_diagram<double>::cell_container_type;
            GeoMapVoronoiDiagramAdaptor(std::vector<GeoMap::GeoNode*>& geo_nodes,std::vector<const GeoMap::GeoNode*>& sites);
            uint64_t get_nodes_in_cell_num(uint64_t cell_index) const;
            std::vector<const GeoMap::GeoNode*> get_nodes_in_cell(uint64_t cell_index) const;
            const std::vector<const GeoMap::GeoNode*> get_nodes_has_osm_tag_in_cell(uint64_t cell_index) const;
            const std::vector<const GeoMap::GeoNode*> get_nodes_has_specified_osm_tag_in_cell(uint64_t cell_index,OSMMapFeature map_feature) const;
            uint64_t get_nodes_has_specified_osm_tag_in_cell_num(uint64_t cell_index,OSMMapFeature map_feature) const;
            std::vector<uint64_t> get_cell_indexes_in_circle_domain(double radius,Point2 loc_utm_xy) const;
            std::vector<uint64_t> get_cell_indexes_which_loc_in(Point2 loc_utm_xy) const;
            void graph_adapt(const voronoi_cell_container_type& cells);
            void shortest_path_distance_to_cells(uint64_t cell_index,std::vector<double>& results) ;
            inline uint64_t sites_num(){ return m_sites_num;}
        private:
            algorithm::Graph m_cell_graph;
            uint64_t m_sites_num;
            std::vector<OSMMapRegion<const GeoMap::GeoNode*>> m_cell_regions;

        };

        }   // namespace map
    } // namespace geo
    
} // namespace name

#endif