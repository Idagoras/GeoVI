#include "algorithm.h"
#include <random>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace geovi::algorithm::distance;
using namespace geovi::algorithm::sample;
using namespace geovi::algorithm::semantic;
using namespace geovi::algorithm::voronoi_diagram;
using namespace geovi::geo::map;
using namespace geovi;

using TagField = enum TagField {
    name = 0,
    map_feature = 1,
    feature_value = 2
};

std::pair<bool,OSMMapFeature> hasFeature(const GeoMap::GeoNode& node,OSMMapFeature feature){
    for(auto feature_tuple : node.features){
        if( get<TagField::map_feature>(feature_tuple) == feature )
            return std::make_pair(true,get<TagField::map_feature>(feature_tuple));
    }
    return std::make_pair(false,OSMMapFeature::None);
}


int64_t loc_num_in_site(int64_t index,std::vector<std::map<std::pair<geo::map::OSMMapFeature,std::string>,int64_t>>& maps){
    auto& map = maps[index];
    int64_t result = 0;
    for(auto it = map.begin(); it != map.end(); ++ it){
        result += it->second;
    }
    return result;
}

// ShortestPathCalculator

void ShortestPathCalculator::calculateShortestPathViaDijkstra(Graph& graph,Index origin_index){
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

// SemanticCategoryCalculator



// PrivacySensitivityCalculator



// CellGrowingAndMergingCalculator

CellGrowingAndMergingCalculator::CellGrowingAndMergingCalculator(
        std::weak_ptr<geovi::algorithm::voronoi_diagram::VoronoiDiagram> vd,
        std::weak_ptr<geovi::geo::map::GeoMap> gmap, double merge_threshold, double grow_threshold):
        m_voronoi_diagram(vd),
        m_gmap(gmap),
        m_merge_threshold(merge_threshold),
        m_grow_threshold(grow_threshold),
        m_psc(gmap)
{
    std::shared_ptr<GeoMap> share_map_ptr = m_gmap.lock();
    std::shared_ptr<VoronoiDiagram> share_vd_ptr = m_voronoi_diagram.lock();
    if( share_map_ptr && share_vd_ptr ){
        m_regions = std::vector<Region>(vd.lock()->sites().capacity());
        auto points = share_map_ptr->getUTMNodesCoordinate();
        for(auto point : points){
            std::pair<Point2,VoronoiDiagram::CellIndex> site_pair = share_vd_ptr->cellIncludePoint(point);
            VoronoiDiagram::CellIndex index = site_pair.second;
            Point2 site = site_pair.first;
            m_regions[index].nodes.push_back(point);
            m_regions[index].value += m_psc.get(point);
            m_regions[index].index = index;
        }
    }
}

struct MergedRegion {
    std::vector<CellGrowingAndMergingCalculator::Region> regions;
    double value;
    void merge(CellGrowingAndMergingCalculator::Region r){
        regions.push_back(r);
        value += r.value;
    }
    CellGrowingAndMergingCalculator::Region merged(){
        CellGrowingAndMergingCalculator::Region r;
        for(auto region : regions){
            merge(r,region);
        }

        return r;
    }
private:
    void merge(CellGrowingAndMergingCalculator::Region& r1,CellGrowingAndMergingCalculator::Region& r2){
        r1.nodes.insert(r1.nodes.end(),r2.nodes.begin(),r2.nodes.end());
        r1.value += r2.value;
    }

};

void regionGrowing(std::vector<MergedRegion>& result,VoronoiDiagram& voronoi_diagram,
                   std::vector<CellGrowingAndMergingCalculator::Region>& regions,
                   double grow_threshold) {

    std::vector<int> visited(regions.capacity(),0);
    for(VoronoiDiagram::CellIndex i = 0; i < regions.capacity(); i++ ){
        auto neighbors_indexes = voronoi_diagram.neighbors(i);
        MergedRegion merge_region ;
        merge_region.merge(regions[i]);
        visited[i] = 1;
        for(VoronoiDiagram::CellIndex k = 0 ; k < merge_region.regions.capacity() ; k ++ ){
            for(VoronoiDiagram::CellIndex neighbor_index : neighbors_indexes){
                if( visited[neighbor_index] == 0 ){
                    double diff = abs(regions[neighbor_index].value-regions[i].value);
                    if( diff < grow_threshold ){
                            merge_region.merge(regions[neighbor_index]);
                            visited[neighbor_index] = 1;
                    }
                }
            }
        }
        result.push_back(merge_region);
    }

}

void regionMerging(std::vector<MergedRegion>& merged_regions,
                   std::vector<CellGrowingAndMergingCalculator::Region>& regions,
                   VoronoiDiagram& voronoi_diagram,
                   double merge_threshold) {

}

