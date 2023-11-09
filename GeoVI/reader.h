#ifndef GEOVI_IO_READER_H
#define GEOVI_IO_READER_H

#include <osmium/io/xml_input.hpp>
#include <osmium/handler.hpp>
#include <string>

namespace geovi
{
    namespace io{
        class Reader {
        public:
            
            osmium::io::Header header;
            Reader(std::string file_name):osmium_reader(file_name,osmium::osm_entity_bits::all){
                header = osmium_reader.header();
                
            }
            osmium::io::Reader& getOSMReader();
           // void readOSMEntityFromHandler(osmium::osm_entity_bits::type osm_ebits,osmium::handler::Hander handler){}
            
        private:
            osmium::io::Reader osmium_reader;
        };

    } // namespace io 
} // namespace GeoVI
#endif