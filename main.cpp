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
#include "ransac.hxx"

void read_SIMPLE(std::vector<shade::point3r>& pts)
{
	shade::point3r x0 = { -1.0, -1.0, -1.0 };
	shade::point3r x1 = { -1.0, -1.0,  1.0 };
	shade::point3r x2 = { -1.0,  1.0, -1.0 };
	shade::point3r x3 = { -1.0,  1.0,  1.0 };
	shade::point3r x4 = { 1.0, -1.0, -1.0 };
	shade::point3r x5 = { 1.0, -1.0,  1.0 };
	shade::point3r x6 = { 1.0,  1.0, -1.0 };
	shade::point3r x7 = { 1.0,  1.0,  1.0 };
	shade::point3r x8 = { 0.0,  0.0,  0.0 };
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

void read_OBJ(const std::string& filename, std::vector<shade::point3r>& pts)
{
	std::ifstream infile(filename);

	std::string line;
	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string type;
		iss >> type;
		if (type == "v")
		{
			shade::point3r pt;
			if (iss >> pt[0] >> pt[1] >> pt[2])
			{
				pts.push_back(pt);
			}
		}
	}
}
void write_to_file(const std::vector<shade::point3r>& pts, std::string filename)
{
	std::ofstream of(filename, std::ofstream::out);
	for (auto it = pts.begin(); it != pts.end(); ++it)
	{
		of << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
	}
	of.close();
}

#include <random>
void create_points(
	const shade::point3r& l,
	const shade::point3r& h,
	const size_t N,
	std::vector<shade::point3r>& pts)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> dx(l[0], h[0]);
	std::uniform_real_distribution<double> dy(l[1], h[1]);
	std::uniform_real_distribution<double> dz(l[2], h[2]);
	pts.reserve(N);
	for (size_t i = 0; i<N; ++i) 
	{
		pts.push_back({ {dx(generator), dy(generator), dz(generator)} });
	}
}

bool test_ztree_canonics()
{
	const shade::point3r l = { {0.0, 0.0, 0.0} };
	const shade::point3r h = { {1.0, 1.0, 1.0} };
	const size_t N = 100000;
	std::vector<shade::point3r> pts;
	create_points(l, h, N, pts);
	using tree_type = shade::ztree<6>;
	tree_type tr(pts);

	const shade::point3r center = { {0.4, 0.4, 0.4} };
	const double radius = 0.4;

	size_t cntr = 0;
	double avg = 0.0;
	double max = 0.0, min = DBL_MAX;
	shade::Sphere sph(center, radius);
	{for (tree_type::iterator it = tr.begin(&sph); it != tr.end(); ++it)
	{
		const double dist = sph.compute_distance_measure(*it);
		avg += dist;
		cntr++;
		max = dist > max ? dist : max;
		min = dist < min ? dist : min;
	}}
	avg /= static_cast<double>(cntr);
	
	std::cout << "Points clashing " << cntr << ", average dist to sphere center " << avg << std::endl;

	return true;
}

template <typename tree_type, typename modl_type>
void test_canonics(
	const tree_type& tr,
	const modl_type& c,
	double& avg,
	double& max,
	double& min,
	size_t& cntr,
	const std::string& directory)
{
	cntr = 0;
	avg = 0.0;
	max = 0.0, min = DBL_MAX;
	std::vector<shade::point3r> pts;
	{for (auto it = tr.begin(&c); it != tr.end(); ++it)
	{
		const double dist = c.compute_distance_measure(*it);
		avg += dist;
		cntr++;
		max = dist > max ? dist : max;
		min = dist < min ? dist : min;
		pts.push_back(*it);
	}}
	avg /= static_cast<double>(cntr);

	//size_t _cntr = 0;
	//double _avg = 0.0;
	//double _max = 0.0, _min = DBL_MAX;
	//tr.get_clashing_pts(&c, pts);
	//{for (auto it = pts.begin(); it != pts.end(); ++it)
	//{
	//	const double dist = c.compute_distance_measure(*it);
	//	_avg += dist;
	//	_cntr++;
	//	_max = dist > _max ? dist : _max;
	//	_min = dist < _min ? dist : _min;
	//}}
	//_avg /= static_cast<double>(_cntr);

	const std::string ofilename(directory + c.get_name() + ".txt");
	write_to_file(pts, ofilename);
}

int main(const int argc, const char * argv[])
{
	//std::vector<shade::point3r> pts;

	//    read_SIMPLE(pts);

	//std::string ifilename("/Users/matthias/Documents/STL/OBJ_Kongming_Wooden_Lock.obj");
	//read_OBJ(ifilename, pts);

	//std::string ofilename("/Users/matthias/Documents/STL/OBJ_Kongming_Wooden_Lock.txt");
	//write_to_file(pts, ofilename);

	//const size_t DEPTH = 5;
	//shade::ztree<DEPTH> t(pts);

	std::vector<shade::point3r> pts;
	{
		const shade::point3r l = { { 0.0, 0.0, 0.0 } };
		const shade::point3r h = { { 1.0, 1.0, 1.0 } };
		const size_t N = 1000000;
		create_points(l, h, N, pts);
	}
	using tree_type = shade::ztree<7>;
	const tree_type tr(pts);

	const shade::point3r center = { { 0.4, 0.4, 0.4 } };
	const shade::point3r normal = { { 0.0, 1.0, 0.0 } };
	const double radius = 0.4;

	//std::shared_ptr<shade::Canonic> sph(new shade::Sphere(center, radius));
	const shade::Sphere sph(center, radius);
	const shade::Plane pln(center, normal);
	const shade::Cylinder cyl(center, normal, radius);
	const shade::Torus tor(center, normal, radius, radius / 4.0);

	std::vector<shade::point3r> sph_pts0, pln_pts0, cyl_pts0, tor_pts0;
	std::vector<shade::point3r> sph_pts1, pln_pts1, cyl_pts1, tor_pts1;

	std::string directory("C:\\Users\\UQ2\\Documents\\Visual Studio 2015\\Projects\\ShaDe\\");
	double avg, max, min;
	size_t cntr;
	//test_canonics(tr, sph, avg, max, min, cntr, directory); // ok
	tr.get_clashing_pts(&sph, sph_pts0);
	tr.get_clashing_pts_greedy(&sph, sph_pts1);

	////test_canonics(tr, pln, avg, max, min, cntr, directory); // ok
	//tr.get_clashing_pts(&pln, pln_pts0);
	//tr.get_clashing_pts_greedy(&pln, pln_pts1);

	////test_canonics(tr, cyl, avg, max, min, cntr, directory); // ok
	//tr.get_clashing_pts(&cyl, cyl_pts0);
	//tr.get_clashing_pts_greedy(&cyl, cyl_pts1);

	////test_canonics(tr, tor, avg, max, min, cntr, directory); // ok
	//tr.get_clashing_pts(&tor, tor_pts0);
	//tr.get_clashing_pts_greedy(&tor, tor_pts1);

	return 0;
}
