#ifndef GEOVI_MECHANISM_H
#define GEOVI_MECHANISM_H

#include <map>
#include <vector>
#include "geomap.h"


namespace geovi {
    namespace algorithm{
        namespace mechanism {
            enum Attack {
                Optimal,
                Beyans
            };

            enum DistanceMeasurement {
                Euclid
            };

            enum ScoreType {
                AE,
                Q_loss,
                PC,
                Time
            };

            class Mechanism {
            public:
                using Attack = enum Attack;
                using DistanceMeasurement = enum DistanceMeasurement;
                using ScoreType = enum ScoreType;
                using Score = std::map<ScoreType,double>;
                using DiscreteDistribution = std::vector<double>;

                virtual void buildDistribution(float _epsilon);
                virtual void computerInferenceFunction();
                virtual void computerAE(Attack attack);
                virtual void computerPC();
                virtual void computerQ_loss();
                Score score();

            protected:
                DiscreteDistribution prior;
                DiscreteDistribution dist;
                float epsilon;
                float start_time;
                float end_time;
            private:
                
            };


            class PLMG : public Mechanism {
            
            };

            class GEM : public Mechanism {
            public:
                GEM(geovi::geo::map::GeoMap& map);

                virtual void buildDistribution(float _epsilon) override;
                virtual void computerInferenceFunction() override;
                virtual void computerAE(Attack attack) override;
                virtual void computerPC() override;
                virtual void computerQ_loss() override;

            private:
                geovi::geo::map::GeoMap& geomap;
            };

            class GVEM : public Mechanism {

            };
        }
    }
}



#endif