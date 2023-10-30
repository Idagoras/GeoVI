#include "GeoVI/convert.h"
#include <iostream>
#include <string>

int main(){
    std::cout << "utem zone = " << std::to_string(int(geovi::LongitudeBands::band_31)) << std::endl;
    geovi::CoordinateSystemConverter csc = geovi::CoordinateSystemConverter(geovi::LongitudeBands::band_31);
    geovi::Point2 point = {49.0,2.0};
    csc.convert(geovi::CoordinateSystemType::WGS84,geovi::CoordinateSystemType::UTM,point);
    std::cout<< "utm zone31 Easting:" << point.x * 10000 << " Northing:" << point.y * 10000 << std::endl;
    return 0;
}