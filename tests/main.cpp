//
//  main.cpp
//  ShaDe
//
//  Created by Matthias Messner on 3/17/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#include <string>

#include <vtkOBJReader.h>
#include <vtkSmartPointer.h>
template <typename T> using vtkptr = vtkSmartPointer<T>;

#include "ztree.hxx"
#include "ransac.hxx"

int main(const int argc, const char * argv[])
{
    const std::string filename("/Users/matthias/Documents/STL/OBJ_Kongming_Wooden_Lock.obj");

    // create reader
    vtkptr<vtkOBJReader> reader = vtkptr<vtkOBJReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // get points
    vtkptr<vtkPoints> points(reader->GetOutput()->GetPoints());
    const vtkIdType npoints = points->GetNumberOfPoints();
 
    // copy points
    std::vector<shade::point_type> pts;
    pts.reserve(npoints);
    for (vtkIdType id=0; id<npoints; ++id) {
        double coords[3];
        points->GetPoint(id, coords);
        pts.push_back({{coords[0], coords[1], coords[2]}});
    }
    
    // build tree
    const size_t DEPTH = 5;
    using tree_type = shade::ztree<DEPTH>;
    tree_type t(pts);
    
    // plane model
    std::array<shade::point_type, shade::Plane::num_pts> rpts;
    t.get_random_points(rpts);
    shade::Plane p(rpts);
    
    // find points on model
    std::vector<tree_type::pt_id_type> pt_ids;
    t.clashing_with(p, pt_ids);
    
    return 0;
}
