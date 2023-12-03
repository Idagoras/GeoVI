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
                public:
                    SemanticCategoryCalculator(geovi::algorithm::voronoi_diagram::VoronoiDiagram& voronoi_diagram,geo::map::GeoMap& g_map);
                    void TF_IDF_UserSemanticsCategory(const UserPOIS& upois,const POIS& pois,geovi::geo::map::GeoMap::GeoNode& loc);
                    void TF_IDF_GeoSemanticsCategory(const geovi::geo::map::GeoMap::GeoNode& node
                                                    ,std::map<geovi::geo::map::OSMMapFeature,double>& features_tf_idf_map
                                                    ,std::map<std::pair<geo::map::OSMMapFeature,std::string>,double>& features_value_tf_idf_map);
                private:
                    std::vector<std::map<geo::map::OSMMapFeature,std::map<std::string,int64_t>>> m_cell_feature_total_occurrences_maps;
                    std::vector<std::map<std::pair<geo::map::OSMMapFeature,std::string>,int64_t>> m_cell_feature_value_total_occurrences_maps;
                    std::map<geo::map::OSMMapFeature,int64_t> m_cell_has_feature_num_map;
                    std::map<std::pair<geo::map::OSMMapFeature,std::string>,int64_t> m_cell_has_feature_value_num_map;
                    geovi::algorithm::voronoi_diagram::VoronoiDiagram& m_voronoi_diagram;
                    geo::map::GeoMap& m_g_map;

                };

                class SemanticSimilarityCalculator {
                    static semantic_similarity calculate(const POI& poi_1,const POI& poi_2);
                };

                class PrivacySensitivityCalculator {
                public:
                    PrivacySensitivityCalculator(std::weak_ptr<geovi::geo::map::GeoMap> gmap);
                    double get(const Point2& p);
                private:
                    std::map<int64_t,double> m_ps_map;
                };

                class CellGrowingAndMergingCalculator{
                public:
                    using Region = struct Region{
                        std::vector<Point2> nodes;
                        double value = 0.0;
                        int64_t index = -1;
                    };
                    CellGrowingAndMergingCalculator(std::weak_ptr<geovi::algorithm::voronoi_diagram::VoronoiDiagram> vd,std::weak_ptr<geovi::geo::map::GeoMap> gmap,double merge_threshold,double grow_threshold);


                private:
                    std::vector<Region> m_regions;
                    std::vector<double> m_semantics;
                    double m_merge_threshold;
                    double m_grow_threshold;
                    PrivacySensitivityCalculator m_psc;
                    std::weak_ptr<geovi::geo::map::GeoMap> m_gmap;
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