#ifndef GEOVI_MECHANISM_H
#define GEOVI_MECHANISM_H

#include <map>
#include <vector>
#include <cmath>
#include "geomap.h"
#include "convert.h"

namespace geovi {
    namespace algorithm{
        namespace mechanism {
            enum Attack {
                Optimal,
                Bayesian
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

                virtual CheckInData computer(const CheckInData& check_in_data);
                virtual Trajectory computer(const Trajectory& trajectory_data);
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
                virtual CheckInData computer(const CheckInData& check_in_data) override;
                virtual void computerInferenceFunction() override;
                virtual void computerAE(Attack attack) override;
                virtual void computerPC() override;
                virtual void computerQ_loss() override;

            private:
                geovi::geo::map::GeoMap& geomap;
            };

            class GVEM : public Mechanism {

            };

            class DP3_SLOC : public Mechanism {
            public:
                DP3_SLOC(double rl,double rm);
                virtual void buildDistribution(float _epsilon) override;
                virtual void computerInferenceFunction() override;
                virtual void computerAE(Attack attack) override;
                virtual void computerPC() override;
                virtual void computerQ_loss() override;
                inline double requirementFunction(double l){ return pow(l/l_s,2); }
                inline double InverseRequirementFunction(double req){ return sqrt(req)*l_s;}

            private:
                std::vector<std::tuple<float,Point2,Point2>> E;
                double r_large;
                double r_small;
                double l_s;

            };
        }
    }
}



#endif