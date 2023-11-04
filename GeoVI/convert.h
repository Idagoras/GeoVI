#ifndef GEOVI_CONVERT_H
#define GEOVI_CONVERT_H

#include <iostream>
#include <proj/coordinateoperation.hpp>
#include <string>

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

using LongtitudeBands= enum LongitudeBands{
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
    Point2(int a,int b):x(a),y(b){}
};

struct Segement{
    Point2 p0;
    Point2 p1;
    Segement(int x1,int y1,int x2,int y2):p0(x1,y1),p1(x2,y2){}
};

struct Line {
    Point2 p0;
    double slope;
    Direction dir = Direction::in;
};

class CoordinateSystemConverter{
public:
    CoordinateSystemConverter(){};
    CoordinateSystemConverter(LongitudeBands bd):band(bd){};
    bool convert(CoordinateSystemType srcCRS,CoordinateSystemType targetCRS,Point2& point);

    
private:
    LongitudeBands band;
    DatabaseContextNNPtr dbContext = DatabaseContext::create();
    bool convertBetweenWGS84AndUTM(Point2& point);
    bool convertBetweenUTMAndWGS84(Point2& point);
    
};
}

#endif