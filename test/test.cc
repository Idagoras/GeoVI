#include "GeoVI/geomap.h"

using namespace geovi::geo::map;
using namespace geovi::io;


int main(){
    OSMReader reader("../OSM/nanjing.xml");
    GeoMap mymap = GeoMap(reader,GeoMap::GeoMapShapeType::circle,GeoMap::Circle());
    return 0;
}