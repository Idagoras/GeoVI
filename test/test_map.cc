//
// Created by 吴伟 on 2023/11/30.
//
#include "GeoVI/geomap.h"
#include "GeoVI/reader.h"
#include "GeoVI/voronoi.h"
#include "GeoVI/algorithm.h"
#include <iostream>
#include <memory>
#include <random>

using namespace std;
using namespace geovi;
using namespace geovi::io;
using namespace geovi::geo::map;
using namespace geovi::algorithm::voronoi_diagram;
using namespace geovi::algorithm::semantic;


void print_node(const GeoMap::GeoNode* geo_node_ptr){
    std::cout << "name is " << geo_node_ptr -> name << std::endl;
    std::cout << "index is " << geo_node_ptr -> index << " osm object id is " << geo_node_ptr -> id << std::endl;
    std::cout << "latitude is " << geo_node_ptr -> loc.latitude << " longitude is " << geo_node_ptr -> loc.longitude << std::endl;
    int k = 0;
    for(auto& feature_tuple : geo_node_ptr -> features){
        std::cout << "map feature " << k << " is " <<  get<0>(feature_tuple) << " and it's value is " << get<2>(feature_tuple) << std::endl;
        ++k;
    }
    std::cout << std::endl;
}

void print_utm_point(const Point2& point){
    std::cout << "point 's utm.x is " << point.x << " and it's utm.y is " << point.y << std::endl;
}

void print_wgs84_point(const Point2& wgs84_point){
    std::cout << "point 's latitude is " << wgs84_point.x << " and it's longitude is " << wgs84_point.y << std::endl;
}

int main(){
    SemanticManager sm;
    std::string bar("bar"),college("college"),ditch("ditch"),anchor("anchor"),bakery("bakery");
    std::cout << "bar-collage : " << sm.semantic_distance(OSMMapFeature::Amenity,bar,OSMMapFeature::Amenity,college) <<std::endl;
    std::cout << "bar-ditch : " << sm.semantic_distance(OSMMapFeature::Amenity,bar,OSMMapFeature::Barrier,ditch) <<std::endl;
    std::cout << "bar-anchor : " << sm.semantic_distance(OSMMapFeature::Amenity,bar,OSMMapFeature::Historic,anchor) <<std::endl;
    std::cout << "bar-bakery : " << sm.semantic_distance(OSMMapFeature::Amenity,bar,OSMMapFeature::Shop,bakery) <<std::endl;
    //sm.iterate();
    /*
    std::shared_ptr<CrossingFilter> filter = make_shared<CrossingFilter>();
    unique_ptr<OSMReader> reader = make_unique<OSMReader>("../OSM/nanjing.xml");
    shared_ptr<GeoMap> geo_map = make_shared<GeoMap>(*reader,filter);

    // 打印顶点数量
    std::cout << "map has " << geo_map->numOfNodes() << " nodes " << std::endl << std::endl;

    // 打印边的数量
    std::cout << "map has " << geo_map->numOfWays() << " ways" << std::endl << std::endl;

    /*
    // 打印所有顶点的WGS84坐标集合
    vector<Point2> wgs84_points = geo_map->getUTMNodesCoordinate();
    int64_t index = 0 ;
    for(auto& wgs84_point : wgs84_points){
        //print_utm_point(wgs84_point);
        ++ index;
    }
    std::cout << std::endl;
*/
    /*
    // 打印地图中所有顶点的信息
    vector<GeoMap::GeoNode*> geo_nodes_ptrs = geo_map->getGeoNodes();
    index = 0;
    for(auto geo_node_ptr : geo_nodes_ptrs){
        std::cout << "point " << index << std::endl;
        print_node(geo_node_ptr);

        auto nodes_ptr = geo_map->find(500,geo_node_ptr->utm_xy.x,geo_node_ptr->utm_xy.y);
        for(auto node : nodes_ptr)
            print_node(node);

        std::cout << std::endl;
        ++ index;
    }
    */
/*
        vector<GeoMap::GeoNode*> geo_nodes_ptrs = geo_map->getGeoNodes();
        auto sites = filter ->crossing_geo_nodes();
        shared_ptr<GeoMapVoronoiDiagramAdaptor> geomap_voronoi_adaptor = make_shared<GeoMapVoronoiDiagramAdaptor>(geo_nodes_ptrs,sites);
        VoronoiDiagram vd;
        VoronoiDiagramBuilder builder(NumericalAccuracy::meter,geo_map->utm_boundary());
        auto sites_utm_points = filter->crossing_utm_xy_points();
        builder.build(vd,sites_utm_points);
        geomap_voronoi_adaptor->graph_adapt(vd.cells());
        std::vector<double> shortest_paths;
        geomap_voronoi_adaptor->shortest_path_distance_to_cells(2,shortest_paths);
        uint64_t cell_index = 0;
        for(auto distance : shortest_paths){
            std::cout << " cell 2 to cell " << cell_index << " shortest path distance = " << distance <<std::endl;
            ++ cell_index;
        }
*/
    /*
        std::default_random_engine e;
        std::uniform_int_distribution<int> u(0,10000);
        std::uniform_int_distribution<int> u_1(0,27);
        std::uniform_int_distribution<int> u_2(0,3);
        e.seed(std::time(0));
        std::vector<GeoMap::GeoNode> nodes;
        std::vector<GeoMap::GeoNode*> nodes_ptr;
        std::vector<const GeoMap::GeoNode*> sites;
        for(int i = 0; i < 100 ; i++ ){
            GeoMap::GeoNode node;
            node.index = i;
            node.utm_xy.x = u(e);
            node.utm_xy.y = u(e);
           // std::cout << "x = " << node.utm_xy.x <<  " y = " << node.utm_xy.y << std::endl;

            for(int j = 0; j <= u_2(e) ; ++ j){
                node.features.push_back({"",OSMMapFeature(u_1(e)),""});
            }
            nodes.push_back(node);


        }

        int index = 0;
        for(auto& node : nodes){
            nodes_ptr.push_back(&node);
            if(index%10 == 0){
                sites.push_back(&node);
            }
            index ++ ;
        }

        for(auto ptr : nodes_ptr){
            std::cout << "x = " << ptr->utm_xy.x <<  " y = " << ptr->utm_xy.y << std::endl;
        }

        shared_ptr<GeoMapVoronoiDiagramAdaptor> geomap_voronoi_adaptor = make_shared<GeoMapVoronoiDiagramAdaptor>(nodes_ptr,sites);
*/
        /*
        int index = 0;
        for(auto distance : distances){
            std::cout << "origin " << origin_ptr->loc.latitude << "," << origin_ptr->loc.longitude << " to " << "vertex " << geo_nodes_ptrs[index]->loc.latitude << ","
            << geo_nodes_ptrs[index]->loc.longitude << " distance is " << distance << std::endl;
            index ++ ;
        }
         */








}