#include "GeoVI/mechanism.h"
#include "GeoVI/algorithm.h"


namespace gam  = geovi::algorithm::mechanism;
namespace ggm  = geovi::geo::map;
namespace gad  = geovi::algorithm::distance;

using namespace geovi::algorithm::sample;


void calculatePriorDistribution(ggm::GeoMap& map,gam::Mechanism::DiscreteDistribution& dist){
    
}

gam::GEM::GEM(geovi::geo::map::GeoMap& map):geomap(map){
    calculatePriorDistribution(map,prior);
}

void gam::GEM::buildDistribution(float _epsilon){
    // gad::ShortestPathCalculator().calculateShortestPath();
    dist = std::vector<double>(prior.capacity());
    int64_t index = 0 ;
    double mass_sum = 0;
    for(auto distance : prior){
        dist[index] = exp(-(_epsilon/2))*distance;
        mass_sum += dist[index];
        index ++ ;
    }
    for(int64_t i = 0; i < index ; i++ ){
        dist[i] /= mass_sum;
    }

}

geovi::CheckInData gam::GEM::computer(const geovi::CheckInData& check_in_data){
    // 扰动到一个维诺区域
    int sample_result = DiscreteDistributionSampler::SampleFromDiscreteDistribution(dist);
    // 从扰动后的维诺区域中选择一个满足语义的点
}
