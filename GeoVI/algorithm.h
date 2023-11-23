#ifndef GEOVI_ALGORITHM_H
#define GEOVI_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "voronoi.h"
#include "geomap.h"
#include "convert.h"


namespace geovi{
    namespace algorithm{

        namespace distance {

            using ShortestPathMatrix = boost::numeric::ublas::matrix<double>;

            class ShortestPathCalculator{
            public:
                using Index = int64_t;
                void calculateShortestPathViaDijkstra(Graph& graph,Index origin_index);
            };
        }

        namespace sample {
            class DiscreteDistributionSampler {
            public:
                using Weights = std::vector<double>;

                static int SampleFromUniformIntDistribution(int min,int max);
                static int SampleFromDiscreteDistribution(Weights weights);
                
            };

            class ContinuousDistributionSampler {

            };
        }

        namespace semantic {
                using role = enum role{
                    user = 0,
                    maintainer = 1
                };
                using stay_time = float;
                using semantic_similarity = float;

                using POI = struct POI {
                    Point2 loc;
                    stay_time st;

                };
                using UserPOIS = std::vector<geovi::geo::map::GeoMap::GeoNode>;
                using POIS = std::vector<geovi::geo::map::GeoMap::GeoNode>;
                class SemanticCategoryCalculator{
                    static void TF_IDF_UserSemanticsCategory(const UserPOIS& upois,const POIS& pois);
                    static void TF_IDF_GeoSemanticsCategory(const POIS& large,const POIS& small);
                };

                class SemanticSimilarityCalculator {
                    static semantic_similarity calculate(const POI& poi_1,const POI& poi_2);
                };

                class PrivacySensitivityCalculator {
                    
                };

                class CellGrowingAndMergingCalculator{
                public:
                    CellGrowingAndMergingCalculator(std::weak_ptr<geovi::algorithm::voronoi_diagram::VoronoiDiagram> vd,double merge_threshold,double grow_threshold):
                    m_voronoi_diagram(vd),
                    m_merge_threshold(merge_threshold),
                    m_grow_threshold(grow_threshold)
                    {};
                    void regionGrowing();
                    void regionMerging();
                private:
                    double m_merge_threshold;
                    double m_grow_threshold;
                    PrivacySensitivityCalculator m_psc;
                    SemanticSimilarityCalculator m_ssc;
                    std::weak_ptr<geovi::algorithm::voronoi_diagram::VoronoiDiagram> m_voronoi_diagram;
                };
        }

        namespace distribution {
                class PlanarLaplacian {

                };

                
        }
    }
}



#endif