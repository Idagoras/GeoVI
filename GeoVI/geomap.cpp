#include "geomap.h"
#include <osmium/handler.hpp>
#include <osmium/osm.hpp>
#include <osmium/io/any_input.hpp>
#include <osmium/visitor.hpp>

class GeoMapHandler : public osmium::handler::Handler{
            public:
                

                GeoMapHandler(geovi::geo::map::GeoMap& map):gmap(map){
                    osmium::handler::Handler();
                };
                void way(const osmium::Way& way);
                void node(const osmium::Node& node);
                void relation(const osmium::Relation& relation);
            private:
                geovi::geo::map::GeoMap& gmap;

            };

// class GeoMap

geovi::geo::map::GeoMap::GeoMap(geovi::io::Reader& reader,GeoMapShapeType type,Shape shape):shape_type(type),mshape(shape){
    GeoMapHandler handler(*this);
    osmium::apply(reader.getOSMReader(),handler);
    reader.getOSMReader().close();
}