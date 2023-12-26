#ifndef GEOVI_WRITER_H
#define GEOVI_WRITER_H
#include <string>
#include <vector>


namespace geovi{
    namespace algorithm{
        namespace semantic{
            struct cluster;
        }
    }
}


namespace geovi{
    namespace io{
        class OSMWriter{
        public:
            void wirte_xml(std::vector<geovi::algorithm::semantic::cluster>& clusters,std::string file_name);
        };
    }
}


#endif