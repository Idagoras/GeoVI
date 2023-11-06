#ifndef GEOVI_ALGORITHM_H
#define GEOVI_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace geovi{
    namespace algorithm{
        namespace distance {

            using ShortestPathMatrix = boost::numeric::ublas:matrix<double>;

            class ShortestPathCaculator{
            public:
                void caculateShortestPath();
            };
        }
    }
}



#endif