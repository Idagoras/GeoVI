#include <iostream>
#include <proj/coordinateoperation.hpp>
#include <proj/crs.hpp>
#include <proj/io.hpp>
#include <proj/util.hpp>
#include <string>

using namespace NS_PROJ::crs;
using namespace NS_PROJ::io;
using namespace NS_PROJ::operation;
using namespace NS_PROJ::util;

namespace geovi{

using LongtitudeBands= enum LongitudeBands{
    band_1 = 1,
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

using degree = float;


struct Point2 {
    double x;
    double y;
};


class CoordinateSystemConverter{
public:
    CoordinateSystemConverter(LongitudeBands lonb,LatitudeBands latb){
        
    }
    

private:
    CRSNNPtr sourceCRS;
    CRSNNPtr targetCRS;
    
};
}