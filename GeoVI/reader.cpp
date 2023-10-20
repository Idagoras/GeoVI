#include "reader.h"

using namespace geovi::io;

osmium::io::Reader& Reader::getOSMReader(){
    return this-> osmium_reader;
}