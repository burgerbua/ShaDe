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
#include <vtkPolyDataMapper.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPointData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>

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
    
    std::set<tree_type::pt_id_type> pt_ids_unique;
    for (auto i : pt_ids) {
        pt_ids_unique.insert(i);
    }
    
    vtkptr<vtkPolyData> polydata = vtkptr<vtkPolyData>::New();
    polydata->SetPoints(points);
    
    unsigned char red[3] = {255, 0, 0};
    unsigned char green[3] = {0, 255, 0};
    vtkptr<vtkUnsignedCharArray> colors = vtkptr<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");
    for (vtkIdType id=0; id<npoints; ++id) {
        if (pt_ids_unique.find(id)==pt_ids_unique.end())
            colors->InsertNextTypedTuple(green);
        else
            colors->InsertNextTypedTuple(red);
    }
    polydata->GetPointData()->SetScalars(colors);

    vtkptr<vtkVertexGlyphFilter> glyphFilter = vtkptr<vtkVertexGlyphFilter>::New();
    glyphFilter->SetInputData(polydata);
    glyphFilter->Update();

    vtkptr<vtkPolyDataMapper> mapper = vtkptr<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyphFilter->GetOutputPort());
//    mapper->SetInputData(polydata);

    vtkptr<vtkActor> actor = vtkptr<vtkActor>::New();
    actor->SetMapper(mapper);
//    actor->GetProperty()->SetPointSize(5);
    
    vtkptr<vtkRenderer> renderer = vtkptr<vtkRenderer>::New();
    vtkptr<vtkRenderWindow> renderWindow = vtkptr<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    vtkptr<vtkRenderWindowInteractor> renderWindowInteractor = vtkptr<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    renderer->AddActor(actor);
//    renderer->SetBackground(.3, .6, .3);
    
    renderWindow->Render();
    renderWindowInteractor->Start();
    
    return 0;
}
