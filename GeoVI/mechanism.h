#ifndef MECHANISM_H
#define MECHANISM_H

#include <map>

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

                virtual void buildDistribution(float _epsilon);
                virtual void computerInferenceFunction();
                virtual void computerAE(Attack attack);
                virtual void computerPC();
                virtual void computerQ_loss();
                Score score();

            protected:
                float epsilon;
                float start_time;
                float end_time;
            private:
                
            };


            class PLMG : public Mechanism {

            };

            class GEM : public Mechanism {

            };

            class GVEM : public Mechanism {

            };
        }
    }
}



#endif