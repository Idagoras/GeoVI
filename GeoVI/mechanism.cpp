#include "GeoVI/mechanism.h"
#include "GeoVI/algorithm.h"
#include <sstream>



namespace gam  = geovi::algorithm::mechanism;
namespace ggm  = geovi::geo::map;
namespace gad  = geovi::algorithm::distance;

using namespace geovi::algorithm::sample;
using namespace gam;
using namespace geovi;

// DataStream

std::string DataStream::next_line() {
    std::string str_line;
    if(!m_fin.eof())
        std::getline(m_fin,str_line);
    return str_line;
}

bool DataStream::eof() {
    return m_fin.eof();
}

bool DataStream::open(){
    if(!m_fin.is_open())
        m_fin.open(m_file_path,std::ios::in);
    return m_fin.is_open();
}

void DataStream::close() {
    m_fin.close();
}

// 特化的 DataParser

bool DataParser<Point2>::can_parse(std::string& string_data) {
    m_parsed_results.clear();
    std::string input_string(string_data);
    std::stringstream ss(input_string);
    std::string word;
    while(ss >> word){
        m_parsed_results.push_back(word);
    }
    if(m_parsed_results.size() == 5){
        m_latest_input_string = input_string;
        return true;
    }else{
        m_parsed_results.clear();
        return false;
    }
}

Point2 DataParser<Point2>::parse(std::string& string_data) {
    if(m_parsed_results.size() == 5 && m_latest_input_string == string_data)
        return Point2{std::stod(m_parsed_results[2]),std::stod(m_parsed_results[3])};
    else{
        if(can_parse(string_data))
            return Point2{std::stod(m_parsed_results[2]),std::stod(m_parsed_results[3])};
        else
            return Point2{0,0};
   }
}

// Result

const std::vector<Point2> &Result::results() const {
    return m_results;
}

const Point2& Result::the_latest_result() const {
    if(!m_results.empty()){
        return m_results.back();
    }
    return m_empty_point;
}

double Result::average_distance_Q_loss() const {
    double q_loss = 0.0;
    uint64_t index = 0;
    for(auto& result : m_results){
        double distance = DistanceCalculator::euclidDistance2D(m_input_locs[index],result);
        q_loss += m_distribution->at(index) * m_prior_distribution->at(index)*distance;
    }
    return q_loss;
}

double Result::average_adversarial_distance_error_AE() const {
    double ae = 0.0;
    u_int64_t index = 0.0;
    for(auto& result : m_results){
        double distance = DistanceCalculator::euclidDistance2D(m_input_locs[index],result);
        ae += m_prior_distribution->at(index) * m_distribution->at(index) * m_adversarial_distribution->at(index) * distance;
    }
    return ae;
}

double Result::performance_criterion_PC() const {
    return (average_distance_Q_loss()+1)/(average_adversarial_distance_error_AE()+1);
}

// Mechanism

Mechanism::Mechanism(){
    m_result.m_prior_distribution = &m_prior;
    m_result.m_distribution = &m_dist;
    m_result.m_adversarial_distribution = &m_adversarial;
}

void Mechanism::pull_data(DataStream &stream, once_finished_call_back call_back) {

}

void Mechanism::commit_result(Point2 result,Point2 input_loc) {
    m_result.m_results.push_back(result);
    m_result.m_input_locs.push_back(input_loc);
}
// 不考虑语义的维诺计算方法

GVEM::GVEM(geovi::geo::map::GeoMapVoronoiDiagramAdaptor &adaptor,float epsilon):m_adaptor(adaptor){
    m_epsilon = epsilon;
}

void GVEM::pull_data(geovi::algorithm::mechanism::DataStream &stream,
                     geovi::algorithm::mechanism::Mechanism::once_finished_call_back call_back) {
    DataParser<Point2> loc_parser;
    if(stream.open()) {
        while (!stream.eof()) {
            if (m_execute_num == 0) {
                m_start_time = TimeUtilHelper::get_current_millis();
                std::string string_line = stream.next_line();
                if(loc_parser.can_parse(string_line)) {
                    Point2 loc = loc_parser.parse(string_line);
                    build_distribution(loc);
                    uint64_t disturbance_cell_index = DiscreteDistributionSampler::SampleFromDiscreteDistribution(
                            m_dist);
                    Point2 result;
                    unsigned long long this_begin_time = TimeUtilHelper::get_current_millis();
                    domain_disturbance(disturbance_cell_index, result, loc);
                    commit_result(result,loc);
                    unsigned long long this_end_time = TimeUtilHelper::get_current_millis();
                    (*call_back)(this_end_time-this_begin_time,m_result);
                }
            }
        }
        m_end_time = TimeUtilHelper::get_current_millis();
        stream.close();
    }
}

void GVEM::build_distribution(geovi::Point2 &loc) {
    auto loc_in_cell_indexes  = m_adaptor.get_cell_indexes_which_loc_in(loc);
    std::vector<double> shortest_paths_distance(m_adaptor.sites_num());
    m_dist.resize(m_adaptor.sites_num(),0);
    uint64_t  cell_index = 0 ;
    if(loc_in_cell_indexes.size() == 1){
        cell_index = loc_in_cell_indexes.front();
    }else if(loc_in_cell_indexes.size() > 1){
        cell_index = loc_in_cell_indexes.front();
    }else
        return;
    m_adaptor.shortest_path_distance_to_cells(cell_index,shortest_paths_distance);
    double total_probability_mass = 0;
    uint64_t i = 0;
    for(auto distance : shortest_paths_distance){
        m_dist[i] = std::exp(-(m_epsilon*distance)/2);
        total_probability_mass += m_dist[i];
        i ++ ;
    }
    i = 0;
    for(auto probability_mass : m_dist){
        m_dist[i] = probability_mass/total_probability_mass;
        i ++ ;
    }
}

void GVEM::domain_disturbance(uint64_t index,Point2& result,Point2& loc) {
    auto nodes = m_adaptor.get_nodes_has_osm_tag_in_cell(index);
    DiscreteDistribution dist(nodes.size(),0);
    uint64_t i = 0;
    double total_probability_mass = 0;
    for(auto node : nodes){
        double distance = DistanceCalculator::euclidDistance2D(loc,node->utm_xy);
        dist[i] = std::exp(-(m_epsilon*distance)/2);
        total_probability_mass += m_dist[i];
        i ++ ;
    }
    i = 0;
    for(auto probability_mass : dist){
        dist[i] = probability_mass/total_probability_mass;
        i ++ ;
    }
    uint64_t result_index = DiscreteDistributionSampler::SampleFromDiscreteDistribution(dist);
    result.x = nodes[result_index]->utm_xy.x;
    result.y = nodes[result_index]->utm_xy.y;
}