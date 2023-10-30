#ifndef GEOVI_VORONOI_H
#define GEOVI_VORONOI_H

#include <boost/polygon/voronoi.hpp>
#include "convert.h"

namespace bp = boost::polygon;

namespace geovi{
    namespace algorithm{
        namespace voronoi_diagram
        {
            class VoronoiDiagram{
            public:
                using VD =  bp::voronoi_diagram<double>;
                friend class VoronoiDiagramBuilder;
                VoronoiDiagram();
            private:
                VD vd;
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