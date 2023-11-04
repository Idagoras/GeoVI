#ifndef MECHANISM_H
#define MECHANISM_H

namespace geovi {
    namespace algorithm{
        namespace mechanism {
            enum Attack {
                Optimal
            };

            enum DistanceMeasurement {
                Euclid
            };

            class Mechanism {
            public:
                using Attack = enum Attack;
                using DistanceMeasurement = enum DistanceMeasurement;
                void buildDistribution(float epsilon);
                void computerInferenceFunction();
                void computerAE(Attack attack);
            private:
            };
        }
    }
}



#endif