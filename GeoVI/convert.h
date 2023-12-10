#ifndef GEOVI_CONVERT_H
#define GEOVI_CONVERT_H

#include <iostream>
#include <proj/coordinateoperation.hpp>
#include <string>
#include <chrono>
#include <cmath>
#include <memory>


using namespace NS_PROJ::crs;
using namespace NS_PROJ::io;
using namespace NS_PROJ::operation;
using namespace NS_PROJ::util;

namespace geovi{


using CoordinateSystemType = enum CoordinateSystemType {
    WGS84,
    UTM,
    Cartesian2,
    Cartesian3
};

using NumericalAccuracy = enum NumericalAccuracy{
    meter = 1,
    decimeter = 10,
    centimeter = 100,
    millimeter = 1000,
};

using LongitudeBands= enum LongitudeBands{
    band_1 = 32601,
    band_2,
    band_3,
    band_4,
    band_5,
    band_6,
    band_7,
    band_8,
    band_9,
    band_10,
    band_11,
    band_12,
    band_13,
    band_14,
    band_15,
    band_16,
    band_17,
    band_18,
    band_19,
    band_20,
    band_21,
    band_22,
    band_23,
    band_24,
    band_25,
    band_26,
    band_27,
    band_28,
    band_29,
    band_30,
    band_31,
    band_32,
    band_33,
    band_34,
    band_35,
    band_36,
    band_37,
    band_38,
    band_39,
    band_40,
    band_41,
    band_42,
    band_43,
    band_44,
    band_45,
    band_46,
    band_47,
    band_48,
    band_49,
    band_50,
    band_51,
    band_52,
    band_53,
    band_54,
    band_55,
    band_56,
    band_57,
    band_58,
    band_59,
    band_60
};

using LatitudeBands = enum LatitudeBands{
    A = 'A',
    B,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    J,
    K,
    L,
    M,
    N,
    O,
    P,
    Q,
    R,
    S,
    T,
    U,
    V,
    W,
    X,
    Y,
    Z
};

using Direction = enum Direction{
    in, 
    out, 
};

using degree = float;


struct Point2 {
    double x;
    double y;
    Point2(){x=0;y=0;};
    Point2(double a,double b):x(a),y(b){}
};

struct Segment{
    Point2 p0;
    Point2 p1;
    Segment(double x1,double y1,double x2,double y2):p0(x1,y1),p1(x2,y2){}
};

struct Line {
    Point2 p0;
    double slope;
    Direction dir = Direction::in;
};

struct CheckInData {
    int64_t user_id;
    std::chrono::system_clock::time_point timePoint;
    double latitude;
    double longitude;
    int64_t location_id;
};

// example:39.984702,116.318417,0,492,39744.1201851852,2008-10-23,02:53:04

struct TrajectoryPoint {
    double latitude;
    double longitude;
    double field;
    double altitude;
    double days_f;
    std::chrono::system_clock::time_point timePoint;
};


using TransportationMode = enum Mv {
    walk,
    bike,
    bus,
    car,
    subway,
    train,
    airplane,
    boat,
    run,
    motorcycle
};

using TransportationModeLabel = struct TransportationModeLabel {
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    TransportationMode mode;
};

struct Trajectory{
    int64_t user_id;
    CoordinateSystemType crs;
    std::vector<TransportationModeLabel> trans_labels;
    std::vector<TrajectoryPoint> tj_points;
};

template<typename T>
struct OSMMapNode{
    typedef T osm_map_node_value_type;
    OSMMapNode<osm_map_node_value_type>* next  = nullptr;
    bool is_link = false;
    OSMMapNode<osm_map_node_value_type>* link;
    osm_map_node_value_type node;

};

template<typename T>
struct OSMTagNode{
    int map_feature;
    uint64_t map_node_num = 0 ;
    typedef T osm_map_node_type;
    OSMTagNode<osm_map_node_type>* next = nullptr;
    std::vector<OSMMapNode<osm_map_node_type>> map_nodes;

};

template<typename T>
struct OSMMapRegion{
    uint64_t  map_node_num = 0;
    typedef T osm_tag_node_type;
    OSMMapRegion<osm_tag_node_type>* next = nullptr;
    OSMMapNode<T>* first_node;
    std::vector<OSMTagNode<osm_tag_node_type>> tag_nodes;
};


class CoordinateSystemConverter{
public:
    CoordinateSystemConverter();
    CoordinateSystemConverter(LongitudeBands bd);
    LongitudeBands utm_identity(double longitude,double latitude);
    void convert(CoordinateSystemType srcCRS,CoordinateSystemType targetCRS,Point2& point);

    
private:
    LongitudeBands band;
    void convertBetweenWGS84AndUTM(Point2& point);
    void convertBetweenUTMAndWGS84(Point2& point);
    
};

class DistanceCalculator{
public:
    inline static double euclidDistance2D(const Point2& p1,const Point2& p2){
        return sqrt(pow((p1.x-p2.x),2)+pow((p1.y-p2.y),2));
    } 
    inline static int double_distance_equal(double distance_1,double distance_2){
        const double epsilon = 1e-4;
        double offset = distance_1-distance_2;
        double offset_abs = std::abs(offset);
        if(offset_abs < epsilon){
            return 0;
        }else{
            if( offset > 0)
                return 1;
            else
                return -1;
        }
    }
    
};

class StringAndTimeConverter{
public:
    static std::chrono::system_clock::time_point convertStringToTime(const std::string& timeStr,const char* fmt);
    static std::string formatTime(const std::chrono::system_clock::time_point& timePoint,const char* fmt);

};

}
#endif