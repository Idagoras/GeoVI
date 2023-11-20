#include "reader.h"

using namespace geovi::io;

osmium::io::Reader& OSMReader::getOSMReader(){
    return this-> osmium_reader;
}