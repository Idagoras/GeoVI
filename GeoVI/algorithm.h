#ifndef GEOVI_ALGORITHM_H
#define GEOVI_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include <functional>
#include "voronoi.h"
#include "geomap.h"
#include "convert.h"
#include <sstream>


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
                    static SemanticManager* getInstance();
                    int semantic_distance(geovi::geo::map::OSMMapFeature feature_1,std::string& feature_1_value,
                                          geovi::geo::map::OSMMapFeature feature_2,std::string& feature_2_value);
                    bool have_semantic_mass(geovi::geo::map::OSMMapFeature feature,std::string feature_value);
                    inline void set_semantic_similar_threshold(int threshold_distance){
                        m_threshold_distance = threshold_distance;
                    }
                private:
                    static SemanticManager* m_instance;
                    SemanticManager();
                    static bool is_load ;
                    std::map<std::string,std::vector<std::string>> m_value_to_path;
                    std::map<std::string ,int> m_distance;
                    int m_threshold_distance = 2;
                };

                struct cluster{
                    cluster(const geo::map::GeoMap::GeoNode * ct){
                        centroid = ct;
                    }
                    void add_element(const geovi::geo::map::GeoMap::GeoNode* element);

                    void add_element(const geovi::geo::map::GeoMap::GeoNode* element,
                                     std::function<bool(std::string& feature_1,std::string& feature_value_1,
                                                        std::string& feature_2,std::string& feature_value_2,
                                                        SemanticManager& sm)> similar_function);


                    int max_radius = 0;
                    int categories_num=0;
                    int size =0;
                    const geo::map::GeoMap::GeoNode * centroid;
                    std::vector<int> elements_num_of_categories;
                    std::vector<std::string> categories;
                    std::vector<const geovi::geo::map::GeoMap::GeoNode*> elements;
                };

                class ClusterCalculator{
                public:
                    ClusterCalculator(int cluster_min_size,int cluster_expected_size,int cluster_categories_num);
                    void calculate(std::vector<cluster>& clusters,
                                   const std::vector<geo::map::GeoMap::GeoNode*>& elements,
                                   const std::vector<const geo::map::GeoMap::GeoNode*>& centroids,
                                   geovi::geo::map::GeoMap& g_map,
                                   std::function<bool(std::string& feature_1,std::string& feature_value_1,
                                                      std::string& feature_2,std::string& feature_value_2,
                                                      SemanticManager& sm)> similar_function);
                private:
                    int m_cluster_min_size;
                    int m_cluster_expected_size;
                    int m_cluster_categories_num;
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