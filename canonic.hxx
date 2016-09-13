//
//  canonic.hxx
//  ShaDe
//
//  Created by Matthias Messner on 3/29/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#ifndef canonic_h
#define canonic_h

#include "ransac.hxx"

namespace shade {

	struct Canonic 
	{
		const double sqrt2 = 1.41421356237;
		const double sqrt3 = 1.73205080757;
		virtual std::string get_name() const = 0;
		virtual bool clashing(const point3r& l, const double width) const = 0;
	};

	class Sphere 
		: public Canonic, public Model<point3r, 4> 
	{
		const point3r center;
		const double radius, radius_sq;
	public:
		static Sphere* make_model(std::vector<point3r>& args)
		{
			
		}

		Sphere(const point3r& _center, const double _radius)
			: center(_center), radius(_radius), radius_sq(_radius*_radius)
		{}
		std::string get_name() const override
		{
			return std::string("Sphere");
		}
		double compute_distance_measure(const point3r& p) const override
		{
			const double cdist = norm(p - center);
			return fabs(radius - cdist);
		}
		bool clashing(const point3r& l, const double width) const override
		{
			double min_dist_sq = 0.0;
			double max_dist_sq = 0.0;
			for (int i = 0; i < 3; ++i) 
			{
				const double ci = center[i];
				const double li = l[i];
				const double hi = li + width;
				if (ci<li) 
				{
					const double min = ci - li;
					min_dist_sq += min * min;
					const double max = ci - hi;
					max_dist_sq += max * max;
				}
				else if (ci>hi) 
				{
					const double min = ci - hi;
					min_dist_sq += min * min;
					const double max = ci - li;
					max_dist_sq += max * max;
				}
			}
			if (radius_sq < min_dist_sq)
			{
				return false;
			}
			if (max_dist_sq > 0.0 && radius_sq > max_dist_sq)
			{
				return false;
			}
			return true;
		}
	};

	class Plane : public Canonic, public Model<point3r, 3>
	{
		const point3r root;
		const point3r normal;
	public:
		Plane(const point3r& _root, const point3r& _normal)
			: root(_root), normal(normalize(_normal))
		{}
		std::string get_name() const override
		{
			return std::string("Plane");
		}
		double compute_distance_measure(const point3r& p) const override
		{
			const point3r vec(p - root);
			return fabs(dot(normal, vec));
		}
		bool clashing(const point3r& l, const double width) const override
		{
			const point3r vecl(l - root);
			const double dotl = dot(normal, vecl);
			if (fabs(dotl) > sqrt3 * width)
			{
				return false;
			}
			const point3r h = { {
				normal[0] * vecl[0] > 0.0 ? 0.0 : width,
				normal[1] * vecl[1] > 0.0 ? 0.0 : width,
				normal[2] * vecl[2] > 0.0 ? 0.0 : width} };
			const point3r vech(l + h - root);
			const double doth = dot(normal, vech);
			return dotl*doth < 0.0;
		}
	};

	class Cylinder : public Canonic, public Model<point3r, 4> 
	{
		const point3r root;
		const point3r axis;
		const double radius;
		point3r get_spine_point(const point3r& p) const
		{
			const point3r vec(p - root);
			const point3r dir(axis * dot(axis, vec));
			return point3r(root + dir);
		}
	public:
		Cylinder(const point3r& _root, const point3r& _axis, const double _radius)
			: root(_root), axis(normalize(_axis)), radius(_radius)
		{}
		std::string get_name() const override
		{
			return std::string("Cylinder");
		}
		double compute_distance_measure(const point3r& p) const override
		{
			const point3r sp = get_spine_point(p);
			const double dist = norm(p - sp);
			return fabs(dist - radius);
		}
		bool clashing(const point3r& l, const double width) const override
		{
			const point3r sp = get_spine_point(l);
			const double dist = norm(l - sp);
			const double diam = sqrt3 * width;
			if (dist + diam < radius || 
				dist - diam > radius) 
			{
				return false;
			}
			const bool inside = dist < radius ? true : false;
			for (int i = 1; i < 8; ++i)
			{
				const point3r h = { {
					l[0] + (i & 1 ? width : 0.0),
					l[1] + (i & 2 ? width : 0.0),
					l[2] + (i & 4 ? width : 0.0)} };
				const point3r sp = get_spine_point(h);
				const double dist = norm(h - sp);
				const bool _inside = dist < radius ? true : false;
				if (inside != _inside)
				{
					return true;
				}
			}
			if (!inside && width>radius*sqrt2)
			{
				return true;
			}
			return false;
		}
	};

	class Torus : public Canonic, public Model<point3r, 4>
	{
		const point3r root;
		const point3r axis;
		const double radius0, radius1;
		point3r get_spine_point(const point3r& p) const
		{
			const point3r vec = cross(point3r(p - root), axis);
			const point3r dir = normalize(cross(axis, vec)) * radius0;
			return point3r(root + dir);
		}
	public:
		Torus(const point3r& _root, const point3r& _axis, const double _radius0, const double _radius1)
			: root(_root), axis(normalize(_axis)), radius0(_radius0), radius1(_radius1)
		{}
		std::string get_name() const override
		{
			return std::string("Torus");
		}
		double compute_distance_measure(const point3r& p) const override
		{
			const point3r sp = get_spine_point(p);
			const double dist = norm(point3r(p - sp));
			return fabs(dist - radius1);
		}
		bool clashing(const point3r& l, const double width) const override
		{
			const point3r sp = get_spine_point(l);
			const double dist = norm(l - sp);
			const double diam = sqrt3 * width;
			if (dist + diam < radius1 ||
				dist - diam > radius1)
			{
				return false;
			}
			const bool inside = dist < radius1 ? true : false;
			for (int i = 1; i < 8; ++i) 
			{
				const point3r h = { {
					l[0] + (i & 1 ? width : 0.0),
					l[1] + (i & 2 ? width : 0.0),
					l[2] + (i & 4 ? width : 0.0)} };
				const point3r sp = get_spine_point(h);
				const double dist = norm(point3r(h - sp));
				const bool _inside = dist < radius1 ? true : false;
				if (inside != _inside)
				{
					return true;
				}
			}
			if (!inside && width>radius1*sqrt2)
			{
				return true;
			}
			return false;
		}
	};

}

#endif /* canonic_h */
