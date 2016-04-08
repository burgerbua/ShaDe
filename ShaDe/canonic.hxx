//
//  canonic.hxx
//  ShaDe
//
//  Created by Matthias Messner on 3/29/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#ifndef canonic_h
#define canonic_h

namespace shade {
    
    struct Canonic {
        virtual bool clashing(const point_type& l, const double width) const = 0;
    };
    
    class Sphere : public Canonic {
        const point_type center;
        const double radius_sq;
    public:
        Sphere(const point_type& _center, const double _radius)
        : center(_center), radius_sq(_radius*_radius) {}
        bool clashing(const point_type& l, const double width) const {
            double min_dist_sq = 0.0;
            double max_dist_sq = 0.0;
            for (int i=0; i<3; ++i) {
                const double ci = center[i];
                const double li = l[i];
                const double hi = li + width;
                if (ci<li) {
                    const double min = ci - li;
                    min_dist_sq += min * min;
                    const double max = ci - hi;
                    max_dist_sq += max * max;
                }
                else if (ci>hi) {
                    const double min = ci - hi;
                    min_dist_sq += min * min;
                    const double max = ci - li;
                    max_dist_sq += max * max;
                }
            }
            if (radius_sq<min_dist_sq) {
                return false;
            }
            if (radius_sq>max_dist_sq) {
                return false;
            }
            return true;
        }
    };
    
    class Plane : public Canonic {
        const point_type root;
        const point_type normal;
    public:
        Plane(const point_type& _root, const point_type& _normal)
        : root(_root), normal(normalize(_normal)) {}
        bool clashing(const point_type& l, const double width) const {
            const point_type vecl(l - root);
            const double dotl = dot(normal, vecl);
            if (fabs(dotl)>3.0*width*width) {
                return false;
            }
            const point_type h = {{
                normal[0]*vecl[0]>0.0 ? 0.0 : width,
                normal[1]*vecl[1]>0.0 ? 0.0 : width,
                normal[2]*vecl[2]>0.0 ? 0.0 : width}};
            const point_type vech(l + h - root);
            const double doth = dot(normal, vech);
            return dotl*doth<0.0;
        }
    };

    class Cylinder : public Canonic {
        const point_type root;
        const point_type axis;
        const double radius;
    public:
        Cylinder(const point_type& _root, const point_type& _axis, const double _radius)
        : root(_root), axis(normalize(_axis)), radius(_radius) {}
        bool clashing(const point_type& l, const double width) const {
            const point_type vecl(l - root);
            const point_type crossl = cross(axis, vecl);
            const double distl = norm(crossl);
            const double diam = sqrt(3.0) * width;
            if (distl+diam < radius || distl-diam > radius) {
                return false;
            }
            const double diffl = distl-radius;
            for (int i=1; i<8; ++i) {
                const point_type h = {{
                    l[0] + (i&1 ? width : 0.0),
                    l[1] + (i&2 ? width : 0.0),
                    l[2] + (i&4 ? width : 0.0)}};
                const point_type vech(h - root);
                const double disth = norm(cross(axis, vech));
                const double diffh = disth-radius;
                if (diffh*diffl<0.0) {
                    return true;
                }
            }
            return false;
        }
    };

    class Torus : public Canonic {
        const point_type root;
        const point_type axis;
        const double radius0, radius1;
    public:
        Torus(const point_type& _root, const point_type& _axis, const double _radius0, const double _radius1)
        : root(_root), axis(normalize(_axis)), radius0(_radius0), radius1(_radius1) {}
        bool clashing(const point_type& l, const double width) const {
            const point_type veclr(l - root);
            const double distlr = norm(veclr);
            const double diam = sqrt(3.0) * width;
            if (distlr+diam < radius0-radius1 || distlr-diam > radius0+radius1) {
                return false;
            }
            const double distln = dot(axis, veclr);
            if (distln+diam < radius1 || distln-diam > radius1) {
                return false;
            }
            const point_type crossl = cross(axis, veclr);
            const point_type dirl = normalize(cross(crossl, axis));
            const point_type rootl = {{
                root[0]+dirl[0]*radius0,
                root[1]+dirl[1]*radius0,
                root[2]+dirl[2]*radius0}};
            const point_type vecl(l - rootl);
            const double distl = norm(vecl);
            const double diffl = distl-radius1;
            for (int i=1; i<8; ++i) {
                const point_type h = {{
                    l[0] + (i&1 ? width : 0.0),
                    l[1] + (i&2 ? width : 0.0),
                    l[2] + (i&4 ? width : 0.0)}};
                const point_type vechr(h - root);
                const point_type crossh = cross(axis, vechr);
                const point_type dirh = normalize(cross(crossh, axis));
                const point_type rooth = {{
                    root[0]+dirh[0]*radius0,
                    root[1]+dirh[1]*radius0,
                    root[2]+dirh[2]*radius0}};
                const point_type vech(l - rooth);
                const double disth = norm(vech);
                const double diffh = disth-radius1;
                if (diffh*diffl<0.0) {
                    return true;
                }
            }
            return true;
        }
    };

}

#endif /* canonic_h */
