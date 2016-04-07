//
//  ztree.hpp
//  ShaDe
//
//  Created by Matthias Messner on 3/21/16.
//  Copyright © 2016 burgerbua. All rights reserved.
//

#ifndef ztree_h
#define ztree_h

#include <array>
#include <cstdint>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "point.hxx"
#include "canonic.hxx"

namespace shade {
    
    uint32_t sep2(uint32_t x) {
        x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
        x = (x ^ (x << 16)) & 0x030000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        return x;
    }
    
    uint32_t cmp2(uint32_t x) {
        x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x >>  8)) & 0x030000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
        return x;
    }
    
    uint32_t enc(uint32_t x, uint32_t y, uint32_t z) {
        return (sep2(x) << 2) | (sep2(y) << 1) | sep2(z);
    }
    
    void dec(uint32_t m, uint32_t& x, uint32_t& y, uint32_t& z) {
        x = cmp2(m >> 2);
        y = cmp2(m >> 1);
        z = cmp2(m);
    }
    
    template <size_t DEPTH>
    class ztree
    {
        const static size_t DIM = 1 << DEPTH;
        
        const std::vector<point_type>& pts;
        
        point_type center;
        double diam;
        
        typedef std::list<size_t> mapped_type;
        typedef std::map<uint32_t, mapped_type::iterator> bst_type;
        mapped_type spts;
        bst_type bst;
        
    public:
        ztree(const std::vector<point_type>& _pts)
        : pts(_pts) {
            const double dbl_max = std::numeric_limits<double>::max();
            point_type lo = { dbl_max,  dbl_max,  dbl_max};
            point_type hi = {-dbl_max, -dbl_max, -dbl_max};
            for (auto it = pts.cbegin(); it != pts.end(); ++it) {
                const point_type& pt = *it;
                if (pt[0]<lo[0]) lo[0] = pt[0];
                if (pt[0]>hi[0]) hi[0] = pt[0];
                center[1] += pt[1];
                if (pt[1]<lo[1]) lo[1] = pt[1];
                if (pt[1]>hi[1]) hi[1] = pt[1];
                center[2] += pt[2];
                if (pt[2]<lo[2]) lo[2] = pt[2];
                if (pt[2]>hi[2]) hi[2] = pt[2];
            }
            diam = 0.0;
            for (int i=0; i<3; ++i) {
                const double curr_diam = hi[i] - lo[i];
                center[i] = lo[i] + curr_diam/2.0;
                if (curr_diam>diam) {
                    diam = curr_diam;
                }
            }
            
            point_type fac;
            for (int i=0; i<3; ++i) {
                fac[i] = center[i] / diam - 0.5;
            }
            
            typedef std::pair<uint32_t, size_t> zpt;
            const size_t npts = pts.size();
            std::vector<zpt> zpts; zpts.reserve(npts);
            for (size_t i = 0; i < npts; ++i) {
                const point_type& pt = pts.at(i);
                const double m0 = pt[0] / diam - fac[0];
                const uint32_t X = static_cast<uint32_t>(std::floor(m0 * (DIM-1)));
                const double m1 = pt[1] / diam - fac[1];
                const uint32_t Y = static_cast<uint32_t>(std::floor(m1 * (DIM-1)));
                const double m2 = pt[2] / diam - fac[2];
                const uint32_t Z = static_cast<uint32_t>(std::floor(m2 * (DIM-1)));
                zpts.push_back(std::make_pair(enc(X, Y, Z), i));
            }
            
            struct zcomp {
                bool operator()(const zpt& l, const zpt& r) const {
                    return l.first < r.first;
                }
            };
            std::sort(zpts.begin(), zpts.end(), zcomp());
            
            uint32_t cell = 0;
            for (auto it = zpts.cbegin(); it != zpts.cend(); ++it) {
                spts.push_back(it->second);
                if (it->first > cell || it == zpts.cbegin()) {
                    cell = it->first;
                    bst[cell] = std::prev(spts.end());
                }
            }
        }
        void clashing_with(const Canonic& c, std::vector<bst_type::iterator>& nodes) const {
            for (auto it: bst) {
                double x, y, z;
                dec(it->first, x, y, z);
                if (c.clashing({{x, y, z}}, diam)) {
                    nodes.push_back(it);
                }
            }
        }
    };
    
} /* namespace shade */

#endif /* ztree_h */
