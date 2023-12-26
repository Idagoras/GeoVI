#ifndef GEOVI_MECHANISM_H
#define GEOVI_MECHANISM_H

#include <map>
#include <vector>
#include <cmath>
#include "geomap.h"
#include "convert.h"
#include <fstream>
#include "algorithm.h"

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

            class DataStream{
            public:
                DataStream(const char* path);
                bool open();
                bool eof();
                void close();
                std::string next_line();
            private:
                std::ifstream m_fin;
                std::string m_file_path;
            };

            template<class T>
            class DataParser{
            public:
                bool can_parse(std::string& string_data);
                T parse(std::string& string_data);
            };

            template<>
            class DataParser<Point2>{
            public:
                bool can_parse(std::string& string_data);
                Point2 parse(std::string& string_data);
            private:

                bool m_is_parsed;
                std::vector<std::string> m_parsed_results;
                std::string m_latest_input_string;
            };

            class Result {
            public:
                friend class Mechanism;
                const std::vector<Point2>& results() const;
                const Point2& the_latest_result() const;
                double average_distance_Q_loss() const;
                double average_adversarial_distance_error_AE() const;
                double performance_criterion_PC() const;
            private:
                const std::vector<double>* m_adversarial_distribution;
                const std::vector<double>* m_distribution;
                const std::vector<double>* m_prior_distribution;
                std::vector<Point2> m_results;
                Point2 m_empty_point = {0,0};
                std::vector<Point2> m_input_locs;
            };


            class Mechanism {
            public:
                using Attack = enum Attack;
                using DistanceMeasurement = enum DistanceMeasurement;
                using ScoreType = enum ScoreType;
                using Score = std::map<ScoreType,double>;
                using DiscreteDistribution = std::vector<double>;
                using once_finished_call_back = void(*)(uint64_t duration_millis,const Result& result);
                Mechanism();
                virtual void pull_data(DataStream& stream,once_finished_call_back call_back);
                inline unsigned long long execute_time_millis(){ return m_end_time - m_start_time; }
                inline double average_execute_time_millis(){ return execute_time_millis()/m_execute_num ;}



            protected:
                void commit_result(Point2 result,Point2 input_loc);
                DiscreteDistribution m_adversarial;
                DiscreteDistribution m_prior;
                DiscreteDistribution m_dist;
                uint64_t m_execute_num = 0;
                float m_epsilon;
                unsigned long long m_start_time;
                unsigned long long m_end_time;
                Result m_result;

            };


            class PLMG : public Mechanism {

            };

            class GEM : public Mechanism {
            public:
                GEM(geovi::geo::map::GeoMap& map);
            private:

            };

            class GVEM : public Mechanism {
            public:
                GVEM(geovi::geo::map::GeoMapVoronoiDiagramAdaptor& adaptor,float epsilon);
                virtual void pull_data(DataStream& stream,once_finished_call_back call_back) override;
            private:
                void build_prior_distribution();
                void build_adversarial_distribution();
                void build_distribution(Point2& loc);
                void domain_disturbance(uint64_t index,Point2& result,Point2& loc);
                geovi::geo::map::GeoMapVoronoiDiagramAdaptor& m_adaptor;
            };

            class DP3_SLOC : public Mechanism {
            public:
                DP3_SLOC(double rl,double rm);
                inline double requirementFunction(double l){ return pow(l/l_s,2); }
                inline double InverseRequirementFunction(double req){ return sqrt(req)*l_s;}

            private:
                std::vector<std::tuple<float,Point2,Point2>> E;
                double r_large;
                double r_small;
                double l_s;

            };

            class GSEM : public Mechanism{
            public:
                GSEM(geovi::geo::map::GeoMap& g_map,geovi::geo::map::MapFeatureFilter& filter);
                void pull_data(DataStream& stream,once_finished_call_back call_back) override;
            private:
                void build_prior_distribution();
                void build_adversarial_distribution();
                void build_distribution(Point2& loc);
                void domain_disturbance(uint64_t index,Point2& result,Point2& loc);

                DiscreteDistribution m_crossing_distribution;
                DiscreteDistribution m_location_distribution;
                geovi::geo::map::GeoMap& m_g_map;
                std::vector<geovi::algorithm::semantic::cluster> m_clusters;
                geovi::geo::map::MapFeatureFilter& m_filter;
            };
        }
    }
}



#endif