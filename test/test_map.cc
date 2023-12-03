//
// Created by 吴伟 on 2023/11/30.
//
#include "GeoVI/geomap.h"
#include "GeoVI/reader.h"
#include <iostream>
#include <memory>

using namespace std;
using namespace geovi;
using namespace geovi::io;
using namespace geovi::geo::map;

int main(){
    unique_ptr<OSMReader> reader = make_unique<OSMReader>("../OSM/nanjing.xml");
    shared_ptr<GeoMap> geo_map = make_shared<GeoMap>(*reader,GeoMap::GeoMapShapeType::circle,GeoMap::Circle());

    // 打印顶点数量
    std::cout << "map has " << geo_map->numOfNodes() << " nodes " << std::endl << std::endl;

    // 打印边的数量
    std::cout << "map has" << geo_map->numOfWays() << "ways" << std::endl << std::endl;
/*
    // 打印所有顶点的WGS84坐标集合
    vector<Point2> wgs84_points = geo_map->getUTMNodesCoordinate();
    int64_t index = 0 ;
    for(auto& wgs84_point : wgs84_points){
        std::cout << "point " << index << " 's latitude is " << wgs84_point.x << " and it's longitude is " << wgs84_point.y << std::endl;
        ++ index;
    }
    std::cout << std::endl;

    // 打印地图中所有顶点的信息
    vector<GeoMap::GeoNode*> geo_nodes_ptrs = geo_map->getGeoNodes();
    index = 0;
    for(auto geo_node_ptr : geo_nodes_ptrs){
        std::cout << "point " << index << std::endl;
        std::cout << "name is " << geo_node_ptr -> name << std::endl;
        std::cout << "index is " << geo_node_ptr -> index << " osm object id is " << geo_node_ptr -> id << std::endl;
        std::cout << "latitude is " << geo_node_ptr -> loc.latitude << " longitude is " << geo_node_ptr -> loc.longitude << std::endl;
        int k = 0;
        for(auto& feature_tuple : geo_node_ptr -> features){
            std::cout << "map feature " << k << " is " <<  get<0>(feature_tuple) << " and it's value is " << get<2>(feature_tuple) << std::endl;
            ++k;
        }
        std::cout << std::endl;
        ++ index;
    }
    */
}