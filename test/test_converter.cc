#include "GeoVI/convert.h"
#include <iostream>
#include <string>

using namespace geovi;

int main(){
    Point2 p{118.7700243,31.9608390};
    CoordinateSystemConverter converter(LongitudeBands::band_50);
    converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,p);
    std::cout << p.x << " " << p.y << std::endl;
}