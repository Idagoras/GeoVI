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

                using semantic_similarity = float;


                using UserPOIS = std::vector<const geovi::geo::map::GeoMap::GeoNode*>;
                using POIS = std::vector<const geovi::geo::map::GeoMap::GeoNode*>;
                class SemanticCategoryCalculator{
                public:
                    static std::vector<geovi::geo::map::OSMMapFeature> semantic_category_frequency_inverse_domain_frequency(geovi::geo::map::GeoMapVoronoiDiagramAdaptor& adaptor,
                                                                                                                            Point2 loc,double r_small,double r_large,
                                                                                                                            std::vector<double>& score_out);

                private:


                };

                class SemanticManager{
                public:
                    SemanticManager();
                    int semantic_distance(geovi::geo::map::OSMMapFeature feature_1,std::string& feature_1_value,
                                          geovi::geo::map::OSMMapFeature feature_2,std::string& feature_2_value);

                private:
                    static bool is_load ;
                    std::map<std::string,std::vector<std::string>> m_value_to_path;
                    std::map<std::string ,int> m_distance;
                };



                class CellGrowingAndMergingCalculator{
                public:
                    using Region = struct Region{
                        std::vector<Point2> nodes;
                        double value = 0.0;
                        int64_t index = -1;
                    };
                    CellGrowingAndMergingCalculator(double merge_threshold,double grow_threshold);


                private:

                };
        }

        namespace distribution {
                class PlanarLaplacian {

                };

                
        }
    }
}



#endif