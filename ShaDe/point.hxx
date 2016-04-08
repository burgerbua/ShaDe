//
//  point.hxx
//  ShaDe
//
//  Created by Matthias Messner on 3/29/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#ifndef point_h
#define point_h

namespace shade {
    
    typedef std::array<double, 3> point_type;
    
    inline point_type operator-(const point_type& l, const point_type& r) {
        return {{l[0]-r[0], l[1]-r[1], l[2]-r[2]}};
    }
    inline point_type operator+(const point_type& l, const point_type& r) {
        return {{l[0]+r[0], l[1]+r[1], l[2]+r[2]}};
    }
    inline double dot(const point_type& l, const point_type& r) {
        return l[0]*r[0] + l[1]*r[1] + l[2]*r[2];
    }
    inline point_type cross(const point_type& l, const point_type& r) {
        return {{l[1]*r[2] - l[2]*r[1], l[2]*r[0] - l[0]*r[2], l[0]*r[1] - l[1]*r[0]}};
    }
    inline double norm(const point_type& v) {
        return sqrt(dot(v, v));
    }
    inline point_type normalize(const point_type& vec) {
        const double len = norm(vec);
        return {{vec[0]/len, vec[1]/len, vec[2]/len}};
    }
}

#endif /* point_h */
