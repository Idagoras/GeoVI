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

SemanticCategoryCalculator::SemanticCategoryCalculator(
        geovi::algorithm::voronoi_diagram::VoronoiDiagram &voronoi_diagram,GeoMap& g_map):
        m_voronoi_diagram(voronoi_diagram),
        m_g_map(g_map){
    m_cell_feature_total_occurrences_maps = std::vector<std::map<geo::map::OSMMapFeature,std::map<std::string,int64_t>>>(voronoi_diagram.sites_num()
                                                                                                                        ,std::map<geo::map::OSMMapFeature,std::map<std::string,int64_t>>());
    m_cell_feature_value_total_occurrences_maps = std::vector<std::map<std::pair<geo::map::OSMMapFeature,std::string>,int64_t>>(voronoi_diagram.sites_num()
            ,std::map<std::pair<geo::map::OSMMapFeature,std::string>,int64_t>());
    CoordinateSystemConverter converter;
    auto nodes = m_g_map.getGeoNodes();
    for( auto& node : nodes ){
        Point2 cd_p{node->loc.longitude,node->loc.latitude};
        converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,cd_p);
        std::pair<Point2,VoronoiDiagram::CellIndex> site_pair = voronoi_diagram.cellIncludePoint(cd_p);
        VoronoiDiagram::CellIndex index = site_pair.second;
        auto cell_feature_total_occurrences_map = m_cell_feature_total_occurrences_maps[index];
        for(auto& feature_tuple : node -> features){
            auto node_feature_pair = std::make_pair(std::get<TagField::map_feature>(feature_tuple),std::get<TagField::feature_value>(feature_tuple));
            if(cell_feature_total_occurrences_map.count(std::get<TagField::map_feature>(feature_tuple)) > 0){
                auto& value_num_map =  cell_feature_total_occurrences_map[std::get<TagField::map_feature>(feature_tuple)];
                if(value_num_map.count(std::get<TagField::feature_value>(feature_tuple)) <= 0){
                    m_cell_has_feature_value_num_map[node_feature_pair] += 1;
                }
                value_num_map[std::get<TagField::feature_value>(feature_tuple)] += 1 ;
            }else{
                m_cell_has_feature_num_map[std::get<TagField::map_feature>(feature_tuple)] += 1;
                std::map<std::string,int64_t> value_num_map;
                value_num_map.insert({std::get<TagField::feature_value>(feature_tuple),1});
                cell_feature_total_occurrences_map.insert({std::get<TagField::map_feature>(feature_tuple),value_num_map});
                m_cell_has_feature_value_num_map[node_feature_pair] += 1;
            }
            m_cell_feature_value_total_occurrences_maps[index][node_feature_pair] += 1;
        }

    }

}

void SemanticCategoryCalculator::TF_IDF_UserSemanticsCategory(const UserPOIS& u_pois,const POIS& pois,GeoMap::GeoNode& loc){
    auto features = loc.features;
    std::map<OSMMapFeature,std::map<std::string,int64_t>> feature_total_occurrences_map;
}

void SemanticCategoryCalculator::TF_IDF_GeoSemanticsCategory(const GeoMap::GeoNode& node,std::map<geovi::geo::map::OSMMapFeature,double>& features_tf_idf_map
        ,std::map<std::pair<geo::map::OSMMapFeature,std::string>,double>& features_value_tf_idf_map){
    CoordinateSystemConverter converter;
    Point2 cd_p{node.loc.longitude,node.loc.latitude};
    converter.convert(CoordinateSystemType::WGS84,CoordinateSystemType::UTM,cd_p);
    auto cell_index_pair = m_voronoi_diagram.cellIncludePoint(cd_p);
    int64_t cell_index = cell_index_pair.second;
    auto& cell_map = m_cell_feature_total_occurrences_maps[cell_index];
    int64_t feature_num = 0; // 总的特征数
    int64_t feature_value_num = 0; // 总的特征值数
    int64_t feature_loc_num = 0; // 具有特征的地点数
    std::map<OSMMapFeature,int64_t> node_feature_num_map ; // 具有节点特征的地点数字典

    for(auto& feature_pair : cell_map){
        ++ feature_num ;
        auto feature_value_num_map = feature_pair.second;
        for(auto& feature_value_pair : feature_value_num_map){
            ++ feature_value_num;
            feature_loc_num += feature_value_pair.second;
            auto feature_exist_pair = hasFeature(node,feature_pair.first);
            if(feature_exist_pair.first){
                node_feature_num_map[feature_exist_pair.second] += 1;
            }
        }
    }


    for(auto& feature_tuple : node.features){
        double node_feature_tf
        = node_feature_num_map[std::get<TagField::map_feature>(feature_tuple)]/m_cell_feature_total_occurrences_maps[cell_index].size();
        double node_feature_idf
        = log(m_voronoi_diagram.sites_num()/m_cell_has_feature_num_map[std::get<TagField::map_feature>(feature_tuple)]);
        double node_feature_tf_idf = node_feature_tf * node_feature_idf;
        features_tf_idf_map[std::get<TagField::map_feature>(feature_tuple)] = node_feature_tf_idf;
        auto node_feature_pair = std::make_pair(std::get<TagField::map_feature>(feature_tuple),std::get<TagField::feature_value>(feature_tuple));
        double node_feature_value_tf = m_cell_feature_value_total_occurrences_maps[cell_index][node_feature_pair]/
                                       loc_num_in_site(cell_index,m_cell_feature_value_total_occurrences_maps);
        double node_feature_value_idf = m_voronoi_diagram.sites_num()/m_cell_has_feature_value_num_map[node_feature_pair];
        double node_feature_value_tf_idf = node_feature_value_tf * node_feature_value_idf;
        features_value_tf_idf_map[node_feature_pair] = node_feature_value_tf_idf;

    }

}

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

