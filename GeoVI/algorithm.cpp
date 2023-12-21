#include "algorithm.h"
#include <random>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <memory>
#include <stack>


#define GEOVI_XML_ATTR_STR "<xmlattr>"


using namespace geovi::algorithm::distance;
using namespace geovi::algorithm::sample;
using namespace geovi::algorithm::semantic;
using namespace geovi::algorithm::voronoi_diagram;
using namespace geovi::geo::map;
using namespace geovi;

static boost::property_tree::ptree s_pt;
static boost::property_tree::ptree s_empty_pt;
bool SemanticManager::is_load = false;

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

std::vector<geovi::geo::map::OSMMapFeature> SemanticCategoryCalculator::semantic_category_frequency_inverse_domain_frequency(geovi::geo::map::GeoMapVoronoiDiagramAdaptor& adaptor,
                                                                                                                             Point2 loc,double r_small,double r_large,
                                                                                                                             std::vector<double>& score_out){
   std::vector<uint64_t> sites_index_in_small_radius_domain = adaptor.get_cell_indexes_in_circle_domain(r_small,loc);
   std::vector<uint64_t> sites_index_in_large_radius_domain = adaptor.get_cell_indexes_in_circle_domain(r_large,loc);
   std::vector<const GeoMap::GeoNode *> nodes_in_small_radius;
   uint64_t small_radius_domain_osm_tag_loc_num = 0;
   std::vector<double> small_radius_domain_osm_tags_scf(29);
   for(auto index : sites_index_in_small_radius_domain){
       auto nodes = adaptor.get_nodes_has_osm_tag_in_cell(index);
       nodes_in_small_radius.insert(nodes_in_small_radius.end(),nodes.begin(),nodes.end());
       small_radius_domain_osm_tag_loc_num += adaptor.get_nodes_in_cell_num(index);
       for(int i =0 ; i< 29 ; i++){
           small_radius_domain_osm_tags_scf[i] += adaptor.get_nodes_has_specified_osm_tag_in_cell_num(index,OSMMapFeature(i));
       }
   }
   for(int i = 0 ; i < small_radius_domain_osm_tags_scf.size() ; i ++ ){
       small_radius_domain_osm_tags_scf[i] /= small_radius_domain_osm_tag_loc_num;
   }

   uint64_t  large_radius_domain_num = sites_index_in_large_radius_domain.size();
   std::vector<double> large_radius_osm_tag_idf(29);
   for(auto index : sites_index_in_large_radius_domain){
       for(int i = 0; i < 29 ; i ++ ){
           int osm_tag_num = adaptor.get_nodes_has_specified_osm_tag_in_cell_num(index,OSMMapFeature(i));
           if(osm_tag_num > 0)
               large_radius_osm_tag_idf[i] += 1;
       }
   }
   for(int i = 0; i< large_radius_osm_tag_idf.size() ; i++){
       large_radius_osm_tag_idf[i] = std::log((large_radius_osm_tag_idf[i]+1)/(large_radius_domain_num+1));
   }

   std::vector<double> small_radius_scf_idf(29);
   for(int i = 0; i < 29 ; i++){
       small_radius_scf_idf[i] = large_radius_osm_tag_idf[i] * small_radius_domain_osm_tags_scf[i];
   }

   std::vector<double> small_radius_osm_tag_node_to_loc_min_distance(29,std::numeric_limits<double>::max());
   for(auto node : nodes_in_small_radius){
       double distance = DistanceCalculator::euclidDistance2D(loc,node->utm_xy);
       for(auto feature_tuple : node->features){
           int feature = get<TagField::map_feature>(feature_tuple);
           small_radius_osm_tag_node_to_loc_min_distance[feature] = DistanceCalculator::double_distance_equal(distance,small_radius_osm_tag_node_to_loc_min_distance[feature]) < 0 ? distance : small_radius_osm_tag_node_to_loc_min_distance[feature];
       }
   }
   score_out.resize(29,0);
   for(int i = 0; i < 29 ; i++){
       score_out[i] = small_radius_scf_idf[i] / (small_radius_osm_tag_node_to_loc_min_distance[i] + 0.1);
   }
   double max_score = std::numeric_limits<double>::min();
   std::vector<OSMMapFeature> features_category;
   for(int i = 0; i<29; i++){
       if(score_out[i] > max_score){
           features_category.clear();
           features_category.push_back(OSMMapFeature(i));
           max_score = score_out[i];
       }else if(score_out[i] == max_score)
           features_category.push_back(OSMMapFeature(i));
   }
    return features_category;
}


// SemanticManager


using namespace boost::property_tree;

void build_tree(const ptree& p_tr,std::vector<std::string>& path,
                std::map<std::string,std::vector<std::string>>& result,
                std::map<geovi::geo::map::OSMMapFeature,int>& distance){
    std::string name;
    int d = std::numeric_limits<int>::min();
    BOOST_FOREACH(const ptree::value_type & node, p_tr){
        if(std::string(node.first.data())== GEOVI_XML_ATTR_STR){
            BOOST_FOREACH(const ptree::value_type & attr ,node.second){
                if(attr.first == "type"){
                    name = std::string(attr.second.data());
                    path.push_back(name);
                }
                if(attr.first == "distance"){
                    std::string distance_str(attr.second.data());
                    if(distance_str == "infinity"){
                        d = std::numeric_limits<int>::max();
                    }else
                        d = std::stoi(std::string(attr.second.data()));
                }
            }

        }
    }
    if( d != std::numeric_limits<int>::min()){
        distance.insert({MapFeatureStringConverter::get(name),d});
    }
    std::cout << std::endl;
    bool is_leaf = true;
    BOOST_FOREACH(const ptree::value_type & node, p_tr){
        if(node.first != GEOVI_XML_ATTR_STR){
            build_tree(node.second,path,result,distance);
            is_leaf = false;
        }
    }
    if(is_leaf){
        result.insert({path[1] + ":" + name,path});
    }
    path.pop_back();
}


SemanticManager::SemanticManager() {

    if(!SemanticManager::is_load){
        try{
            read_xml("../GeoVI/semantic_tree.xml",s_pt);
            SemanticManager::is_load=true;
            std::vector<std::string> path;
            path.push_back("root");
            build_tree(s_pt.get_child("semantic-tree.root",s_pt),
                       path,m_value_to_path,m_distance);

            for(auto& pair : m_value_to_path){
                std::cout << pair.first << " : " ;
                for(auto& p : pair.second){
                    std::cout << p << " . " ;
                }
                std::cout << std::endl;
            }
        }catch(const ptree_error& e){
            std::cout << "Error: " << e.what() << std::endl;
        }
    }
}

int SemanticManager::semantic_distance(geovi::geo::map::OSMMapFeature feature_1, std::string &feature_1_value,
                                       geovi::geo::map::OSMMapFeature feature_2, std::string &feature_2_value) {

}




// CellGrowingAndMergingCalculator



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

void regionGrowing(std::vector<MergedRegion>& result,
                   GeoMapVoronoiDiagramAdaptor& adaptor,
                   std::vector<CellGrowingAndMergingCalculator::Region>& regions,
                   double grow_threshold) {

    std::vector<int> visited(regions.capacity(),0);



}


