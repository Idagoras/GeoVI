#ifndef GEOVI_VORONOI_H
#define GEOVI_VORONOI_H

#include <boost/polygon/voronoi.hpp>
#include <boost/numeric/ublas/matrix.hpp>
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
                using matrix_element_type = double;
                using coordinate_type = double;
                using VD =  bp::voronoi_diagram<coordinate_type>;
                using Points = std::vector<Point2>;
                using Segments = std::vector<Segment>;
                using Lines = std::vector<Line>;
                using Matrix = boost::numeric::ublas::matrix<matrix_element_type>;
                using CellIndex = int64_t;
                friend class VoronoiDiagramBuilder;
                VoronoiDiagram(){};
                inline int64_t sites_num(){return m_points.capacity();};
                Points vertices();
                Points sites();
                std::vector<CellIndex> neighbors(CellIndex index);
                Segments finiteEdges();
                Lines infiniteEdges();
                std::pair<Point2,CellIndex> cellIncludePoint(Point2 p);
                Matrix shortestPathBetweenCells();

                
                
                
            private:
                VD m_vd;
                Points m_points;
                Matrix m_shortest_paths;
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

        } // namespace voronoi_diagram
        
    } // namespace algorithm
}

#endif