#ifndef GEOVI_VORONOI_H
#define GEOVI_VORONOI_H

#include <boost/polygon/voronoi.hpp>
#include "convert.h"
#include <vector>

namespace bp = boost::polygon;

template <>
struct bp::geometry_concept<geovi::Point2> { typedef point_concept type; };
  
template <>
struct bp::point_traits<geovi::Point2> {
        typedef int coordinate_type;
   
        static inline coordinate_type get(const geovi::Point2& point, orientation_2d orient) {
            return (orient == HORIZONTAL) ? point.x : point.y;
        }
};      

namespace geovi{
    namespace algorithm{
        namespace voronoi_diagram
        {

        
            class VoronoiDiagram{
            public:
                using VD =  bp::voronoi_diagram<double>;
                using Points = std::vector<Point2>;
                friend class VoronoiDiagramBuilder;
                VoronoiDiagram(){};
                Points vertices();

            private:
                VD vd;
                Points points;
            };

            class VoronoiDiagramBuilder{
            public :
                using Points = std::vector<Point2>;
                using Origin = Point2;
                using InputCoordinateSystemType = geovi::CoordinateSystemType;
                VoronoiDiagramBuilder(NumericalAccuracy numericalAccuracy,Origin originPoint);
                void build(InputCoordinateSystemType inputType,VoronoiDiagram& voronoi_dragram,Points& points);
            private :
                int multiple = 0;  
                Origin origin;
                CoordinateSystemConverter converter;                
            };

        } // namespace voronoi_dragram
        
    } // namesapce algorithm
}

#endif