#ifndef GEOVI_GEOMAP_H
#define GEOVI_GEOMAP_H

#include <osmium/handler.hpp>
#include <osmium/osm.hpp>
#include "reader.h"

namespace geovi
{
    namespace geo
    {
        namespace map{
            class GeoMap{
            public:
                typedef enum {
                    circle,
                    Rectangle,
                } GeoMapShapeType;
                
                typedef struct{
                    // longitude
                        double longitude;
                    // latitude
                        double latitude;

                } Location;

                typedef struct Shape
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
                

                GeoMapShapeType shape_type;
                GeoMap(geovi::io::Reader reader,GeoMapShapeType type,Shape shape);


            private:
                

            };
            class GeoMapHandler : public osmium::handler::Handler{
            public:
                void way(const osmium::Way& way);
                void node(const osmium::Node& node);
                void relation(const osmium::Relation& relation);

            };

            class GeoVoronoiMap{

            };
            
        }   // namespace map
    } // namespace geo
    
} // namespace name

#endif