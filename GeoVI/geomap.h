#ifndef GEOVI_GEOMAP_H
#define GEOVI_GEOMAP_H

#include "reader.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/polygon/voronoi.hpp>
#include "algorithm.h"
#include <vector>



namespace geovi
{
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
                    Millitary, // This is used for facilities and on land used by the military.
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
                
                typedef struct {
                    std::string name;
                    map_object_id_type id;
                    double semantic_sensitivity = -1;
                    OSMMapFeature map_feature = OSMMapFeature::None;
                } GeoNode;

                typedef struct {
                    GeoNode& source;
                    GeoNode& target;
                    map_object_id_type id;
                    std::string name;
                    double capacity;
                } GeoWay;

                GeoMapShapeType shape_type;
                Shape mshape;
                GeoMap(geovi::io::Reader& reader,GeoMapShapeType type,Shape shape);

                bool addNode(GeoNode node);
                bool addWay(GeoNode source,GeoNode target);

                std::vector<Point2> getNodes();
                



               // GeoMap POISOfAnAreaWithCenter(Location loc);
               // GeoMap POISOfAnAreaWithCenterObjectID(map_object_id_type id);


            private:
                geovi::algorithm::Graph graph;
                int nodes_num;
                int ways_num;
                int relations_num;
                float max_lat;
                float max_lon;
                float min_lat;
                float min_lon;
            };   

            

            class GeoVoronoiMap{
            public:
                using VoronoiDiagram = boost::polygon::voronoi_diagram<double>;
            private:
                VoronoiDiagram vd;
                
            };
            
        }   // namespace map
    } // namespace geo
    
} // namespace name

#endif