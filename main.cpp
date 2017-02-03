//
//  main.cpp
//  ShaDe
//
//  Created by Matthias Messner on 3/17/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#include <iostream>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

#include "ztree.hxx"

void read_SIMPLE(std::vector<shade::point_type>& pts) {
    shade::point_type x0 = {-1.0, -1.0, -1.0};
    shade::point_type x1 = {-1.0, -1.0,  1.0};
    shade::point_type x2 = {-1.0,  1.0, -1.0};
    shade::point_type x3 = {-1.0,  1.0,  1.0};
    shade::point_type x4 = { 1.0, -1.0, -1.0};
    shade::point_type x5 = { 1.0, -1.0,  1.0};
    shade::point_type x6 = { 1.0,  1.0, -1.0};
    shade::point_type x7 = { 1.0,  1.0,  1.0};
    shade::point_type x8 = { 0.0,  0.0,  0.0};
    pts.push_back(x0);
    pts.push_back(x1);
    pts.push_back(x2);
    pts.push_back(x3);
    pts.push_back(x4);
    pts.push_back(x5);
    pts.push_back(x6);
    pts.push_back(x7);
    pts.push_back(x8);
}

void read_OBJ(const std::string& filename, std::vector<shade::point_type>& pts) {
    std::ifstream infile(filename);
    
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "v") {
            shade::point_type pt;
            if (iss >> pt[0] >> pt[1] >> pt[2]) {
                pts.push_back(pt);
            }
        }
    }
}
void write_to_file(const std::vector<shade::point_type>& pts, std::string filename) {
    std::ofstream of(filename, std::ofstream::out);
    for (auto it = pts.begin(); it!=pts.end(); ++it) {
        of << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
    }
    of.close();
}

int main(const int argc, const char * argv[])
{
    std::vector<shade::point_type> pts;
    
//    read_SIMPLE(pts);
    
    std::string ifilename("/Users/matthias/Documents/STL/OBJ_Kongming_Wooden_Lock.obj");
    read_OBJ(ifilename, pts);

    std::string ofilename("/Users/matthias/Documents/STL/OBJ_Kongming_Wooden_Lock.txt");
    write_to_file(pts, ofilename);

    const size_t DEPTH = 5;
    shade::ztree<DEPTH> t(pts);
    
    return 0;
}
