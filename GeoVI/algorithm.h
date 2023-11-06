#ifndef GEOVI_ALGORITHM_H
#define GEOVI_ALGORITHM_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <vector>

using namespace boost;

namespace boost{
    enum vertex_semantic_sensitivity_t{
        vertex_semantic_sensitivity = 2042
    };
    BOOST_INSTALL_PROPERTY(vertex,semantic_sensitivity);

    enum vertex_location_t{
        vertex_location = 1111
    };
    BOOST_INSTALL_PROPERTY(vertex,location);

}


namespace geovi{
    namespace algorithm{

        using VertexDescriptor = adjacency_list_traits<vecS,vecS,directedS>::vertex_descriptor;
        using EdgeDescriptor = adjacency_list_traits<vecS,vecS,directedS>::edge_descriptor;


        using VertexProperties = property<vertex_name_t,std::string,
                        property<vertex_index_t,int64_t,
                        property<vertex_semantic_sensitivity_t,double,
                        property<vertex_location_t,int,
                        property<vertex_predecessor_t, VertexDescriptor,
                        property<vertex_distance_t, double>>>>>>;

        using EdgeProperties =  property<edge_name_t,std::string,
                                property<edge_index_t,int64_t,
                                property<edge_weight_t,double>>>;

        using Graph = adjacency_list<vecS,vecS,directedS,VertexProperties,EdgeProperties>;



        using VertexNameMap = property_map<Graph,vertex_name_t>::type ;
        using ConstVertexNameMap = property_map<Graph,vertex_name_t>::const_type ;  
        using VertexIndexMap = property_map<Graph,vertex_index_t>::type ;
        using ConstVertexIndexMap = property_map<Graph,vertex_index_t>::const_type ;
        using VertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::type;
        using ConstVertexSemanticSensitivityMap = property_map<Graph,vertex_semantic_sensitivity_t>::const_type;
        using VertexLocationMap = property_map<Graph,vertex_location_t>::type;
        using ConstVertexLocationMap = property_map<Graph,vertex_location_t>::const_type;


        using EdgeNameMap = property_map<Graph,edge_name_t>::type ;
        using ConstEdgeNameMap = property_map<Graph,edge_name_t>::const_type ;
        using EdgeIndexMap = property_map<Graph,edge_index_t>::type;
        using ConstEdgeIndexMap = property_map<Graph,edge_index_t>::const_type;
        using EdgeCapacityMap = property_map<Graph,edge_capacity_t>::type;
        using ConstEdgeCapacityMap = property_map<Graph,edge_capacity_t>::const_type;


        namespace distance {

            using ShortestPathMatrix = boost::numeric::ublas:matrix<double>;

            class ShortestPathCaculator{
            public:
                void caculateShortestPath();
            };
        }

        namespace sample {
            class DiscreteDistributionSampler {
            public:
                using Weights = std::vector<double>;

                DiscreteDistributionSampler(){};
                int SampleFromUniformIntDistribution(int min,int max);
                int SampleFromDiscreteDistribution(Weights weights);
                
            };

            class ContinuousDistributionSampler {

            };
        }
    }
}



#endif