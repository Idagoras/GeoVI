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
                using CellIndex = uint64_t;

                friend class VoronoiDiagramBuilder;

                VoronoiDiagram()= default;;
                inline uint64_t sites_num(){return m_points.size();};
                Points vertices() const;
                Points sites() const;
                std::vector<CellIndex> neighbors(CellIndex index) const;
                Segments finiteEdges() const;
                Lines infiniteEdges() const;
                std::pair<Point2,CellIndex> cellIncludePoint(Point2 p) const;
                std::vector<Point2> cellPolygon(CellIndex index) const;
                Matrix shortestPathBetweenCells() const;
                std::vector<CellIndex> cells_include_or_on_its_edge(Point2& p) const;


                
                
                
            private:
                VD m_vd;
                Points m_points;
                Matrix m_shortest_paths;

                std::vector<Point2> m_boundary;
            };

            class VoronoiDiagramBuilder{
            public :
                using Points = std::vector<Point2>;
                using Origin = Point2;
                VoronoiDiagramBuilder(NumericalAccuracy numericalAccuracy,const std::vector<Point2>& rect_boundary,Origin originPoint = {0,0});
                void build(VoronoiDiagram& voronoi_diagram,Points& points);
            private :
                int m_multiple = 0;
                Origin m_origin;
                double m_min_x,m_max_x,m_min_y,m_max_y;
            };

        } // namespace voronoi_diagram
        
    } // namespace algorithm
}

#endif