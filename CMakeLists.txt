cmake_minimum_required(VERSION 3.2)
project(GeoVI)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(Osmium REQUIRED COMPONENTS xml)
if(NOT OSMIUM_FOUND)
    message(WARNING "Libosmium not found!\n")
endif()

find_package(EXPAT REQUIRED)
if(NOT EXPAT_FOUND)
    message(WARNING "Libexpat not found!\n")
endif()


include_directories(
    ${OSMIUM_INCLUDE_DIRS}
    ${EXPAT_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/include
)
link_directories(${CMAKE_SOURCE_DIR}/lib)