//
//  main.cpp
//  ShaDe
//
//  Created by Matthias Messner on 3/17/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#include <string>

#include <chrono>

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

using namespace shade;

int main(const int argc, const char * argv[])
{
//    const std::string filename("E:\\dev\\ShaDe\\files\\OBJ_Golfing.obj");
    const std::string filename("/Users/matthias/Code/ShaDe/tests/models/D_396.obj");
    
    // create reader
    vtkptr<vtkOBJReader> reader = vtkptr<vtkOBJReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // get points
    vtkptr<vtkPoints> points(reader->GetOutput()->GetPoints());
    vtkIdType npoints = points->GetNumberOfPoints();
    printf("num pts %u\n", npoints);
    
    const double base_area = 1.0;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> disr(0, 1);
    
    vtkptr<vtkCellArray> polys(reader->GetOutput()->GetPolys());
    polys->InitTraversal();
    vtkptr<vtkIdList> ids = vtkptr<vtkIdList>::New();
    while (polys->GetNextCell(ids)) {
        assert(ids->GetNumberOfIds()==3);
        point_type vtxs[3];
        for (int i=0; i<ids->GetNumberOfIds(); ++i)
            points->GetPoint(ids->GetId(i), vtxs[i].data());
        const shade::point_type vec0(vtxs[1]-vtxs[0]);
        const shade::point_type vec1(vtxs[2]-vtxs[0]);
        const double area = norm(cross(vec0, vec1));
//        printf("%f\n", area);
        if (area>base_area) {
            const int num = static_cast<int>(area/base_area);
            for (size_t i=0; i<num; ++i) {
                const double r[2] = {disr(gen), disr(gen)};
                const double a = 1.0-sqrt(r[0]);
                const double b = sqrt(r[0]) * (1.0-r[1]);
                const double c = sqrt(r[0]) * r[1];
                const shade::point_type A(vtxs[0] * a);
                const shade::point_type B(vtxs[1] * b);
                const shade::point_type C(vtxs[2] * c);
                const shade::point_type P = A + B + C;
                points->InsertNextPoint(P.data());
            }
        }
    }
    
    npoints = points->GetNumberOfPoints();
    printf("num pts %u\n", npoints);
    
    // copy points
    std::vector<shade::point_type> pts;
    pts.reserve(npoints);
    for (vtkIdType id=0; id<npoints; ++id) {
        double coords[3];
        points->GetPoint(id, coords);
        pts.push_back({{coords[0], coords[1], coords[2]}});
    }
    
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	// build tree
    const size_t DEPTH = 6;
    using tree_type = shade::ztree<DEPTH>;
    tree_type t(pts);
    
    // plane model
    std::array<shade::point_type, shade::Plane::num_pts> rpts;
    t.get_random_points(rpts);
    shade::Plane p(rpts);
    
    // find points on model
    std::vector<size_t> pt_ids;
    t.clashing_with(p, pt_ids);
    
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    std::set<size_t> pt_ids_unique;
    for (auto i : pt_ids) {
        pt_ids_unique.insert(i);
    }
	assert(pt_ids.size() == pt_ids_unique.size());

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
