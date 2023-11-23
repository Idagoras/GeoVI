#ifndef GEOVI_EXPERIMENT_H
#define GEOVI_EXPERIMENT_H

#include "mechanism.h"
#include "reader.h"
#include <stdarg.h>
#include <string>


namespace geovi
{
    template<typename T>
    class DataSet {
    public:
        using DataList = std::vector<T>;
        using iterator = typename DataList::iterator;
        using const_iterator = typename DataList::const_iterator;
        const_iterator begin() { return m_ds.begin();}
        const_iterator end() { return m_ds.end(); }
        void description();
    private:
        std::string m_name;
        std::string m_publisher;
        std::vector<T> m_ds;
        int size;


    };

    template<typename T>
    class Data{
    public:
        T real_obj;
        Data(const char* dataStr){ };
    };

    template<typename T>
    class DataSetMaker{
    public:
        using LoadingMethod = enum LoadingMethod {
            local_file,
            download
        };
        using FileType = enum FileType {
            txt,
            rdf
        };
        DataSetMaker(LoadingMethod lm,const char* url,FileType file_type,Data<T> dataType);
        void make(DataSet<T>& ds);
    
    };


    template<typename T>
    class Experiment {
    public:
        using ExperimentResult = enum ExperimentResult {
            AE = 1,
            Q_loss = 2,
            PC = 4,
        };
        
        explicit Experiment(geovi::algorithm::mechanism::Mechanism baseline,geovi::algorithm::mechanism::Mechanism main,DataSet<T> dataSet);
        void start();

    };

    
} // namespace geovi




#endif