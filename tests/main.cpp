//
//  main.cpp
//  ShaDe
//
//  Created by Matthias Messner on 3/17/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#include <string>

#include "ztree.hxx"

#include <vtkOBJReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkSmartPointer.h>
template <typename T> using vtkptr = vtkSmartPointer<T>;

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
    shade::ztree<DEPTH> t(pts);
    
    return 0;
}
