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
}