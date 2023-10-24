cmake_minimum_required(VERSION 3.2)
project(GeoVI)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# libosmium

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(Osmium REQUIRED COMPONENTS xml)
if(NOT OSMIUM_FOUND)
    message(WARNING "Libosmium not found!\n")
endif()

# expat

#find_package(EXPAT REQUIRED)
#if(NOT EXPAT_FOUND)
#    message(WARNING "Libexpat not found!\n")
#endif()

# boost

find_package(Boost REQUIRED)
if(NOT Boost_FOUND)
    message(WARNING "Boost not found!\n")
endif()
message(STATUS "${Boost_LIBRARIES}")

include_directories(
    ${Boost_INCLUDE_DIRS}
    ${OSMIUM_INCLUDE_DIRS}
    ${EXPAT_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

aux_source_directory( GeoVI SRC_LIST)
add_library(geovi ${SRC_LIST})
link_libraries(
    ${BOOST_LIBRARYDIR}
)


add_executable(Test test/test.cc)
target_link_libraries(Test geovi  ${OSMIUM_XML_LIBRARIES})