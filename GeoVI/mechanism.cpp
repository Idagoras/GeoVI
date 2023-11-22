#include "GeoVI/mechanism.h"
#include "GeoVI/algorithm.h"

namespace gam  = geovi::algorithm::mechanism;
namespace ggm  = geovi::geo::map;
namespace gad  = geovi::algorithm::distance;



void caculatePriorDistribution(ggm::GeoMap& map,gam::Mechanism::DiscreteDistribution& dist){
    
}

gam::GEM::GEM(geovi::geo::map::GeoMap& map):geomap(map){
    caculatePriorDistribution(map,prior);
}

void gam::GEM::buildDistribution(float _epsilon){
    // gad::ShortestPathCaculator().caculateShortestPath();
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

}
