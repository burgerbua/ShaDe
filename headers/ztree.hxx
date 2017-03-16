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
#include <set>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "point.hxx"
#include "canonic.hxx"

namespace shade {
    
    inline uint32_t sep2(uint32_t x) {
        x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
        x = (x ^ (x << 16)) & 0x030000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x <<  8)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x <<  4)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x <<  2)) & 0x09249249; // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        return x;
    }
    
    inline uint32_t cmp2(uint32_t x) {
        x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        x = (x ^ (x >>  2)) & 0x030c30c3; // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x >>  4)) & 0x0300f00f; // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x >>  8)) & 0x030000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
        return x;
    }
    
    inline uint32_t enc(uint32_t x, uint32_t y, uint32_t z) {
        return (sep2(x) << 2) | (sep2(y) << 1) | sep2(z);
    }
    
    inline void dec(uint32_t m, uint32_t& x, uint32_t& y, uint32_t& z) {
        x = cmp2(m >> 2);
        y = cmp2(m >> 1);
        z = cmp2(m);
    }
    
    template <size_t DEPTH>
    class ztree
    {
        const static uint32_t MASK = (1 << DEPTH) - 1;
        
        const std::vector<point_type>& pts;
        
        point_type center;
        double width;
        
        typedef std::list<size_t> mapped_type;
        typedef std::map<uint32_t, mapped_type::const_iterator> bst_type;
        mapped_type spts;
        bst_type bst;
        
        void clashing_with(const Canonic& c, const bst_type::const_iterator nl, const bst_type::const_iterator nh, const size_t level, std::vector<bst_type::const_iterator>& nodes) const {
            uint32_t T0, T1, T2;
            dec(nl->first, T0, T1, T2);
            const uint32_t mask = 0xffffffff<<(DEPTH-level);
            T0 &= mask;
            T1 &= mask;
            T2 &= mask;
            const double t0 = static_cast<double>(T0) / static_cast<double>(MASK);
            const double t1 = static_cast<double>(T1) / static_cast<double>(MASK);
            const double t2 = static_cast<double>(T2) / static_cast<double>(MASK);
            const point_type l = {{
                width*t0 + center[0] - width/2.0,
                width*t1 + center[1] - width/2.0,
                width*t2 + center[2] - width/2.0}};
            const double curr_width = width / (1<<level);
            if (c.clashing(l, curr_width)) {
                if (level==DEPTH) {
                    nodes.push_back(nl);
                }
                else {
                    const size_t clevel = level + 1;
                    const uint32_t C = 1<<(DEPTH-clevel);
                    for (short o=0; o<8; ++o) {
                        const uint32_t T0l = o&1 ? T0|C : T0;
                        const uint32_t T1l = o&2 ? T1|C : T1;
                        const uint32_t T2l = o&4 ? T2|C : T2;
                        const uint32_t Ml = enc(T0l, T1l, T2l);
                        auto nlc = std::lower_bound(nl, nh, Ml, [](const bst_type::value_type& n, const uint32_t& _Ml){return n.first<_Ml;});
                        const uint32_t T0h = T0l|C;
                        const uint32_t T1h = T1l|C;
                        const uint32_t T2h = T2l|C;
                        const uint32_t Mh = enc(T0h, T1h, T2h);
                        auto nhc = std::upper_bound(nl, nh, Mh, [](const uint32_t& _Mh, const bst_type::value_type& n){return _Mh<n.first;});
                        clashing_with(c, nlc, nhc, clevel, nodes);
                    }
                }
            }
        }
//        void clashing_with_greedy(const Canonic& c, std::vector<bst_type::iterator>& nodes) const {
//            for (auto it: bst) {
//                double x, y, z;
//                dec(it->first, x, y, z);
//                if (c.clashing({{x, y, z}}, width)) {
//                    nodes.push_back(it);
//                }
//            }
//        }
//        void clashing_with(const Canonic& c, std::vector<bst_type::iterator>& nodes) const {
//            const size_t level = 0;
//            clashing_with(c, bst.begin(), bst.end(), level, nodes);
//        }
        
        uint32_t g2l(const point_type& pt) const
        {
            const double t0 = (pt[0] - center[0]) / width + 0.5;
            const uint32_t T0 = static_cast<uint32_t>(std::floor(t0 * MASK));
            const double t1 = (pt[1] - center[1]) / width + 0.5;
            const uint32_t T1 = static_cast<uint32_t>(std::floor(t1 * MASK));
            const double t2 = (pt[2] - center[2]) / width + 0.5;
            const uint32_t T2 = static_cast<uint32_t>(std::floor(t2 * MASK));
            const uint32_t m = enc(T0, T1, T2);
            return m;
        }
        point_type l2g(const uint32_t m) const
        {
            uint32_t T0, T1, T2;
            dec(m, T0, T1, T2);
            const double t0 = static_cast<double>(T0) / static_cast<double>(MASK);
            const double t1 = static_cast<double>(T1) / static_cast<double>(MASK);
            const double t2 = static_cast<double>(T2) / static_cast<double>(MASK);
            const point_type g = {{
                width*t0 + center[0] - width/2.0,
                width*t1 + center[1] - width/2.0,
                width*t2 + center[2] - width/2.0}};
            return g;
        }
        
    public:
        //using pt_id_type = mapped_type::value_type;
        
        ztree(const std::vector<point_type>& _pts)
        : pts(_pts) {
            const double dbl_max = std::numeric_limits<double>::max();
            point_type lo = {{ dbl_max,  dbl_max,  dbl_max}};
            point_type hi = {{-dbl_max, -dbl_max, -dbl_max}};
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
            width = 0.0;
            for (int i=0; i<3; ++i) {
                const double curr_diam = hi[i] - lo[i];
                center[i] = lo[i] + curr_diam/2.0;
                if (curr_diam>width) {
                    width = curr_diam;
                }
            }
            
            point_type fac;
            for (int i=0; i<3; ++i) {
                fac[i] = center[i] / width - 0.5;
            }
            
            typedef std::pair<uint32_t, size_t> zpt;
            const size_t npts = pts.size();
            std::vector<zpt> zpts; zpts.reserve(npts);
            
            size_t idx = 0;
            for (auto pt : pts) {
                const double t0 = pt[0] / width - fac[0];
                const uint32_t T0 = static_cast<uint32_t>(std::floor(t0 * MASK));
                const double t1 = pt[1] / width - fac[1];
                const uint32_t T1 = static_cast<uint32_t>(std::floor(t1 * MASK));
                const double t2 = pt[2] / width - fac[2];
                const uint32_t T2 = static_cast<uint32_t>(std::floor(t2 * MASK));
                zpts.emplace_back(enc(T0, T1, T2), idx++);
            }
            
            std::sort(zpts.begin(), zpts.end(), [](const zpt& l, const zpt& r){return l.first < r.first;});
            
            uint32_t cell = 0;
            for (auto it = zpts.cbegin(); it != zpts.cend(); ++it) {
                spts.push_back(it->second);
                if (it->first > cell || it == zpts.cbegin()) {
                    cell = it->first;
                    bst[cell] = std::prev(spts.end());
                }
            }
        }

        template <size_t num_pts>
        void get_random_points(std::array<point_type, num_pts>& rpts) const
        {
            // get random leaf node
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<size_t> disn(0, pts.size()-1);
            const size_t id = disn(gen);
            const uint32_t rm = g2l(pts.at(id));
            auto cn = std::prev(std::lower_bound(std::begin(bst), std::end(bst), rm, [](const bst_type::value_type& n, const uint32_t& _rm){return n.first<_rm;}));

            // point counter
            size_t cntr = 0;
            
            // find random points in random leaf
            std::set<size_t> ids;
            auto nn = std::next(cn);
            auto end = (nn==std::end(bst)) ? std::end(spts) : nn->second;
            const size_t size = std::distance(cn->second, end);
            std::uniform_int_distribution<size_t> disp(0, size-1);
            while (ids.size()<size) {
                // use points only once
                const size_t id = disp(gen);
                if (ids.insert(id).second == false)
                    continue;
                auto pt = cn->second;
                std::advance(pt, id);
                rpts[cntr++] = pts.at(*pt);
                // return if enough points found
                if (cntr == num_pts)
                    return;
            }

            // not enough points found, retry
            get_random_points(rpts);
        }
        
        //void clashing_with(const Canonic& c, std::vector<mapped_type::value_type>& pt_ids) const {
        //    // find all nodes clashing with c
        //    std::vector<bst_type::const_iterator> nodes;
        //    const size_t level = 0;
        //    bst_type::const_iterator nl = bst.begin();
        //    bst_type::const_iterator nh = bst.end();
        //    clashing_with(c, nl, nh, level, nodes);
        //    // get all point ids from clashing nodes
        //    for (auto cn : nodes) {
        //        auto nn = std::next(cn);
        //        auto end = (nn==std::end(bst)) ? std::end(spts) : nn->second;
        //        for (auto it = cn->second; it != end; ++it) {
        //            pt_ids.push_back(*it);
        //        }
        //    }
        //}

		void clashing_with(const Canonic& c, std::vector<size_t>& pt_ids) const {
			std::stack<std::pair<uint32_t, int>> todo;
			todo.emplace(0, 0);
			while (!todo.empty()) {
				auto p = todo.top();
				todo.pop();
				const int level = p.second;
				const uint32_t level_mask = 0xffffffff << 3 * (DEPTH - level);
				const uint32_t m = p.first & level_mask;
				const double curr_width = width / (1 << level);
				const point_type l = l2g(m);
				if (c.clashing(l, curr_width)) {
					if (level == DEPTH) {
						auto cn = bst.find(m);
						if (cn != std::end(bst)) {
							auto nn = std::next(cn);
							auto end = (nn == std::end(bst)) ? std::end(spts) : nn->second;
							for (auto id = cn->second; id != end; ++id)
								pt_ids.push_back(*id);
						}
					}
					else {
						const int clevel = level + 1;
						for (int o = 0; o < 8; ++o) {
							const uint32_t cm = m | (o << 3 * (DEPTH - clevel));
							todo.emplace(cm, clevel);
						}
					}
				}
			}
		}


    };
    
} /* namespace shade */

#endif /* ztree_h */
