#include "algorithm.h"
#include <random>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <memory>
#include <stack>
#include <bitset>
#include <algorithm>


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

void get_path(std::string& path_str,std::vector<std::string>& path){
    int index = 0;
    path_str.clear();
    for(auto p : path){
        path_str.append(p);
        if( index != path.size() - 1)
            path_str.append(":");
    }
}

void build_tree(const ptree& p_tr,std::vector<std::string>& path,
                std::map<std::string,std::vector<std::string>>& result,
                std::map<std::string,int>& distance){
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
                        try{
                            d = std::stoi(std::string(attr.second.data()));
                        }catch(std::out_of_range& e_1){
                            std::cout << "Error : distance string : " << attr.second.data() << " convert to int fail >> " << e_1.what();
                            d = std::numeric_limits<int>::min();
                        }catch(std::invalid_argument& ia){
                            std::cout << "Error : distance string : " << attr.second.data() << " convert to int fail >> " << ia.what();
                            d = std::numeric_limits<int>::min();
                        }
                }
            }

        }
    }
    std::string path_str;
    get_path(path_str,path);
    if( d != std::numeric_limits<int>::min()){
        distance.insert({path_str,d});
    }
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

            for(auto pair : m_distance){
                std::cout << " key: " << pair.first << ":distance = " << pair.second << std::endl;
            }
        }catch(const ptree_error& e){
            std::cout << "Error: " << e.what() << std::endl;
        }
    }
}

int SemanticManager::semantic_distance(geovi::geo::map::OSMMapFeature feature_1, std::string &feature_1_value,
                                       geovi::geo::map::OSMMapFeature feature_2, std::string &feature_2_value) {
    std::string key_1 = MapFeatureStringConverter::get(feature_1) + ":" +  feature_1_value;
    std::string key_2 = MapFeatureStringConverter::get(feature_2) + ":" + feature_2_value;
    std::string key_distance_1 = "root:" + MapFeatureStringConverter::get(feature_1) + ":";
    std::string key_distance_2 = "root:" + MapFeatureStringConverter::get(feature_2) + ":";
    int distance_1 = 0,distance_2 = 0;
    try{
        distance_1 = m_distance.at(key_distance_1);
        distance_2 = m_distance.at(key_distance_2);
        auto& path_1 = m_value_to_path.at(key_1);
        auto& path_2 = m_value_to_path.at(key_2);
        std::string path_str_1,path_str_2;
        get_path(path_str_1,path_1);
        get_path(path_str_2,path_2);
        if( distance_1 == std::numeric_limits<int>::max()){
            if(m_distance.find(path_str_1) != m_distance.end())
                distance_1 = m_distance.at(path_str_1);
        }
        if( distance_2 == std::numeric_limits<int>::max()){
            if(m_distance.find(path_str_2) != m_distance.end())
                distance_2 = m_distance.at(path_str_2);
        }
        if(distance_1 == std::numeric_limits<int>::max() || distance_2 == std::numeric_limits<int>::max())
            return std::numeric_limits<int>::max();

        int min_len = std::min(path_1.size(),path_2.size());
        int common_prefix_len = 0;
        while(common_prefix_len < min_len && path_1[common_prefix_len] == path_2[common_prefix_len]){
            common_prefix_len ++ ;
        }
        if(common_prefix_len == 0)
            return std::numeric_limits<int>::max();
        return distance_1 * (path_1.size()-common_prefix_len) + distance_2 * (path_2.size() - common_prefix_len);
    }catch(std::out_of_range& e){
        std::cout  << "Error : " << "key 1 = "<<key_1 << " key 2 = " << key_2 << "key distance 1 = "<<key_distance_1 << " key distance 2 = " << key_distance_2 << " >> " << e.what() << std::endl;
        return -1;
    }
}


// ClusterCalculator

ClusterCalculator::ClusterCalculator(int cluster_min_size, int cluster_expected_size, int cluster_categories_num):m_cluster_min_size(cluster_min_size),
m_cluster_expected_size(cluster_expected_size),m_cluster_categories_num(cluster_categories_num){}


bool check_cluster_satisfy_condition(int min_size,int categories_num,int expected_size,cluster& cl){
    //static const double variance_threshold = 3;

    if( cl.categories_num < categories_num)
        return false;
    if(cl.size >= expected_size)
        return true;
    /*
    double variance = StatisticUtilHelper::variance(cl.elements_num_of_categories);
    if(variance < variance_threshold)
        return true;
    */
    return false;
}


bool check_element_can_be_added_to_the_set(const GeoMap::GeoNode* element,cluster& cl,int cluster_categories_num,int cluster_min_size, int cluster_expected_size){
    if( element->features.empty())
        return false;
    if(cl.size < cluster_min_size)
        return true;
    int max_category_num = std::numeric_limits<int>::min();
    std::vector<std::string> max_num_categories;
    int index = 0;
    for(auto& d : cl.elements_num_of_categories){
        if( d > max_category_num){
            max_category_num = d;
            max_num_categories.clear();
            max_num_categories.push_back(cl.categories[index]);
        }
        if( d == max_category_num )
            max_num_categories.push_back(cl.categories[index]);
        index ++;
    }


    bool new_category = true;
    std::string element_category = std::get<TagField::name>(element->features[0])+":"+std::get<TagField::feature_value>(element->features[0]);
    for(auto& category : cl.categories){
     //   std::cout << "point feature : " << element_category << " cluster category : " << category << std::endl;
        if( category == element_category){
            new_category = false;
            break;
        }
    }
    int cost_current,cost_other;
    if(new_category){
        if(cl.categories_num + 1 <= cluster_categories_num)
            return true;
        else{
            cost_current = max_category_num;
            cost_other = cluster_categories_num;
        }
    }else{
        cost_current = cl.categories_num;
        cost_other = max_category_num;
    }

    if( cost_current < cost_other )
        return true;
    return false;
}


void ClusterCalculator::calculate(std::vector<cluster> &clusters,
                                  const std::vector<geo::map::GeoMap::GeoNode *> &elements,
                                  const std::vector<const geo::map::GeoMap::GeoNode*>& centroids,
                                  geovi::geo::map::GeoMap& g_map,
                                  std::function<bool(const geo::map::GeoMap::GeoNode*,const geo::map::GeoMap::GeoNode*)> similar_function) {
    std::vector<u_int8_t> elements_visited(elements.size(),0);
    std::map<int64_t,u_int8_t> centroids_visited;
    for(auto& centroid : centroids){
        centroids_visited.insert({centroid->index,0});
    }
    clusters.clear();
    int64_t index = 0;
    for(auto& centroid : centroids){
        if(!centroids_visited[centroid->index]){
            std::cout << "centroid " << index << " :" << std::endl;
            clusters.emplace_back(centroid);
            cluster& cl = clusters.back();
            bool completed = false;
            int multi = 1;
            int radius = 100;
            while(!completed){
                auto point_around = g_map.find(radius*multi,centroid->utm_xy.x,centroid->utm_xy.y);
                std::cout << "around radius : " << multi * radius << " has " << point_around.size() << " nodes " << std::endl;
                std::sort(point_around.begin(),point_around.end(),[centroid](const GeoMap::GeoNode* v1,const GeoMap::GeoNode* v2){
                    double distance_1 = DistanceCalculator::euclidDistance2D(v1->utm_xy,centroid->utm_xy);
                    double distance_2 = DistanceCalculator::euclidDistance2D(v2->utm_xy,centroid->utm_xy);
                    return DistanceCalculator::double_distance_equal(distance_1,distance_2) < 0;
                });
                for(auto& point : point_around){
                    if(elements_visited[point->index] == 0){
                        if(check_element_can_be_added_to_the_set(point,cl,m_cluster_categories_num,m_cluster_min_size,m_cluster_expected_size)){
                            cl.add_element(point);
                            elements_visited[point->index] = 1;
                            if( centroids_visited.find(point->index) != centroids_visited.end())
                                centroids_visited[point->index] = 1;
                            completed = check_cluster_satisfy_condition(m_cluster_min_size,m_cluster_categories_num,m_cluster_expected_size,cl);
                            if(completed)
                                break;
                        }
                    }
                }
                completed = completed  || point_around.size() == elements.size();
                if(completed)
                    cl.max_radius = multi*radius;
                multi ++ ;
            }
            index ++;
        }
    }
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


