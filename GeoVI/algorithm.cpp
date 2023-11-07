#include "algorithm.h"
#include <random>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace geovi::algorithm::distance;
using namespace geovi::algorithm::sample;

// ShortestPathCaculator

void ShortestPathCaculator::caculateShortestPathViaDijkstra(Graph& graph,Index origin_index){
    VertexDescriptor origin_vertex = vertex(origin_index,graph);
    VertexPredecessorMap p = get(vertex_predecessor,graph);
    VertexDistanceMap d = get(vertex_distance,graph);
    dijkstra_shortest_paths(graph,origin_vertex,predecessor_map(p).distance_map(d));
}

// Sample

// DiscreteDistributionSampler

int DiscreteDistributionSampler::SampleFromUniformIntDistribution(int min,int max){
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> distribute(min,max);
    return distribute(gen);
}

int DiscreteDistributionSampler::SampleFromDiscreteDistribution(Weights weights){
    std::mt19937 gen(std::random_device{}());
    std::discrete_distribution<int> dist(weights.begin(),weights.end());
    return dist(gen);
}