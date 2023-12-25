#include "writer.h"
#include <iostream>
#include "algorithm.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
using namespace geovi::io;

void OSMWriter::wirte_xml(std::vector<geovi::algorithm::semantic::cluster>& clusters, std::string file_name) {
    // 创建一个空的 PropertyTree
    boost::property_tree::ptree pt;

    // 添加根节点
    boost::property_tree::ptree& rootNode = pt.add("root", "");

    for(auto& cl:clusters){
        boost::property_tree::ptree cl_tree;
        for(auto& element : cl.elements){
            boost::property_tree::ptree el_tree;
            el_tree.put("latitude",element->loc.latitude);
            el_tree.put("longitude",element->loc.longitude);
            cl_tree.add_child("location",el_tree);
        }
        rootNode.add_child("cluster",cl_tree);
    }


    // 写入到 XML 文件
    boost::property_tree::write_xml(file_name, pt);

    std::cout << "XML file created successfully." << std::endl;
}