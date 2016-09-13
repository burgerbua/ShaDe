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
#include <queue>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>

#include "point.hxx"
#include "canonic.hxx"

#define DEBUG

namespace shade {

	inline uint32_t sep2(uint32_t x)
	{
		x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
		x = (x ^ (x << 16)) & 0x030000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x << 8)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x << 4)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x << 2)) & 0x09249249;  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		return x;
	}

	inline uint32_t cmp2(uint32_t x)
	{
		x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
		x = (x ^ (x >> 2)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
		x = (x ^ (x >> 4)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
		x = (x ^ (x >> 8)) & 0x030000ff;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
		x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
		return x;
	}

	inline uint32_t enc(const point3i& lpt)
	{
		return (sep2(lpt[0]) << 2) | (sep2(lpt[1]) << 1) | sep2(lpt[2]);
	}

	inline point3i dec(uint32_t m) 
	{
		uint32_t x = cmp2(m >> 2);
		uint32_t y = cmp2(m >> 1);
		uint32_t z = cmp2(m);
		return {{x, y, z}};
	}


	template <size_t DEPTH>
	class ztree
	{
		const static uint32_t MASK = (1 << DEPTH);

		const std::vector<point3r>& pts;

		double leaf_diam;

		using index_list = std::list<size_t>;
		using node_map = std::map<uint32_t, index_list::iterator>;
		index_list ids;
		node_map bst;

		class Mapper
		{
			point3r lo;
			double diam;

		public:
			double init(
				const point3r& _lo,
				const point3r& _hi)
			{
				// compute center and max diam
				diam = 0.0;
				point3r center;
				for (int i = 0; i < 3; ++i)
				{
					const double curr_diam = _hi[i] - _lo[i] + 2.0 * std::numeric_limits<double>::epsilon();
					center[i] = _lo[i] + curr_diam / 2.0;
					if (curr_diam > diam)
					{
						diam = curr_diam;
					}
				}

				// set low based on center and max diam
				for (int i = 0; i < 3; ++i)
				{
					lo[i] = center[i] - diam / 2.0;
				}
				return diam;
			}
			double get_diam(const size_t level) const
			{
				return diam / static_cast<double>(1 << level);
			}
			point3i glo2loc(
				const point3r& gpt) const
			{
				point3i lpt;
				lpt[0] = static_cast<uint32_t>(std::floor(((gpt[0] - lo[0]) / diam) * MASK));
				lpt[1] = static_cast<uint32_t>(std::floor(((gpt[1] - lo[1]) / diam) * MASK));
				lpt[2] = static_cast<uint32_t>(std::floor(((gpt[2] - lo[2]) / diam) * MASK));
				return lpt;
			}
			point3r loc2glo(
				const point3i& lpt) const
			{
				point3r gpt;
				gpt[0] = lo[0] + diam * (static_cast<double>(lpt[0]) / static_cast<double>(MASK));
				gpt[1] = lo[1] + diam * (static_cast<double>(lpt[1]) / static_cast<double>(MASK));
				gpt[2] = lo[2] + diam * (static_cast<double>(lpt[2]) / static_cast<double>(MASK));
				return gpt;
			}
		} mapper;

		void init_mapper()
		{
			// get low and high of AABB of points
			const double dbl_max = std::numeric_limits<double>::max();
			point3r lo = { dbl_max,  dbl_max,  dbl_max };
			point3r hi = { -dbl_max, -dbl_max, -dbl_max };
			for (auto it = pts.cbegin(); it != pts.end(); ++it)
			{
				const point3r& pt = *it;
				if (pt[0] < lo[0]) lo[0] = pt[0];
				if (pt[0] > hi[0]) hi[0] = pt[0];
				if (pt[1] < lo[1]) lo[1] = pt[1];
				if (pt[1] > hi[1]) hi[1] = pt[1];
				if (pt[2] < lo[2]) lo[2] = pt[2];
				if (pt[2] > hi[2]) hi[2] = pt[2];
			}

			// set mapping object
			leaf_diam = mapper.init(lo, hi) / static_cast<double>(MASK);
		}

		void init()
		{
			init_mapper();

			// compute Z-index for all points
			using zpt = std::pair<uint32_t, size_t>;
			const size_t npts = pts.size();
			std::vector<zpt> zpts; zpts.reserve(npts);
			for (size_t i = 0; i < npts; ++i)
			{
				const point3i lpt = mapper.glo2loc(pts.at(i));
				zpts.push_back(std::make_pair(enc(lpt), i));
			}

			// sort points based on Z-index
			struct zcomp
			{
				bool operator()(const zpt& l, const zpt& r) const
				{
					return l.first < r.first;
				}
			};
			std::sort(zpts.begin(), zpts.end(), zcomp());

			// fill BSP with points based on Z-indices of points
			// all points with same Z-index are in same leaf
			uint32_t cell = 0;
			for (auto it = zpts.cbegin(); it != zpts.cend(); ++it)
			{
				ids.push_back(it->second);
				if (it->first > cell || it == zpts.cbegin())
				{
					cell = it->first;
					bst[cell] = std::prev(ids.end());
				}
			}
		}

	public:
		ztree(const std::vector<point3r>& _pts)
			: pts(_pts)
		{
			init();
		}

		class iterator : public std::iterator<std::forward_iterator_tag, point3r>
		{
			const ztree *const tree;
			const Canonic* c;

			using id_iter = ztree::index_list::const_iterator;
			using nd_iter = ztree::node_map::const_iterator;

			id_iter iit, end_iit;
			nd_iter nit, end_nit;

			bool is_at_node_end() const
			{
				assert(end_iit != iit);
				assert(end_nit != nit);
				nd_iter next_nit = nit;
				if (end_nit == ++next_nit)
				{
					return false;
				}
				return next_nit->second == iit;
			}
			bool is_node_clashing() const
			{
				assert(c);
				const point3i lpt = dec(nit->first);
				if (c->clashing(tree->mapper.loc2glo(lpt), tree->leaf_diam))
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			void advance_if_necessary()
			{
				if (!is_at_node_end())
				{
					return;
				}
				while (end_nit!=nit && (!is_node_clashing() || is_at_node_end()))
				{
					iit = (++nit == end_nit) ? end_iit : nit->second;
				}
			}

		public:
			iterator(const ztree* _tree)
				: tree(_tree)
			{
				assert(!tree->bst.empty());
				end_nit = tree->bst.end();
				nit = end_nit;
				end_iit = tree->ids.end();
				iit = end_iit;
			}
			void init(const Canonic* _c)
			{
				c = _c;

				if (tree->bst.empty())
				{
					return;
				}

				nit = tree->bst.begin();
				iit = nit->second;
				while (end_nit != nit && !is_node_clashing())
				{
					iit = (++nit == end_nit) ? end_iit : nit->second;
				}
			}
			iterator& operator++()
			{
				iit++;
				if (end_iit != iit)
				{
					advance_if_necessary();
				}
				return *this;
			}
			iterator operator++(int) 
			{
				iterator tmp(*this); 
				operator++(); 
				return tmp; 
			}
			bool operator==(const iterator& other) 
			{
				return iit == other.iit; 
			}
			bool operator!=(const iterator& other) 
			{
				return !(*this == other);
			}
			const point3r& operator*() const
			{
				return static_cast<const point3r&>(tree->pts.at(*iit)); 
			}
		};

		iterator begin(const Canonic* c) const 
		{
			iterator it(this);
			it.init(c);
			return it;
		}
		iterator end() const
		{
			return iterator(this);
		}

		void get_clashing_pts(
			const Canonic* c,
			std::vector<point3r>& cpts) const
		{
			cpts.clear();

			struct node_info
			{
				node_map::const_iterator it;
				const size_t level;
				node_info(node_map::const_iterator _it, size_t _level) : it(_it), level(_level) {}
			};

			auto iend = ids.cend();
			auto nend = bst.cend();

			std::queue<node_info> todo;
			todo.emplace(node_info(bst.begin(), 0));
			while (!todo.empty())
			{
				// get last item and pop
				const node_info info = todo.front();
				todo.pop();
				assert(info.it != nend);

				// get lower corner of current box at given level
				const uint32_t mask = 0xffffffff << 3 * (DEPTH - info.level);
				const uint32_t M = info.it->first & mask;
				const point3i lpt = dec(M);

				// if not clashing there is nothing to do
				const double curr_diam = mapper.get_diam(info.level);
				if (info.level > 2 && false == c->clashing(mapper.loc2glo(lpt), curr_diam)) // level>2 only because sphere clashing is bad
				{
					continue;
				}

				// handle clashing nodes
				if (DEPTH == info.level)
				{
					// leaf, hence get all points
					auto nit = info.it;
					auto ciend = (++nit != nend) ? nit->second : iend;
					{for (auto it = info.it->second; it != ciend; ++it)
					{
						cpts.push_back(pts.at(*it));
					}}
				}
				else
				{
					// get all child nodes and add to todo queue
					const size_t clevel = info.level + 1;
					const size_t shift = 3 * (DEPTH - clevel);
					for (short o = 0; o < 8; ++o) {
						const uint32_t _cM = M + (o << shift);
						const uint32_t cM = M | (o << shift);
						auto lo = bst.lower_bound(cM);
						if (lo == nend)
						{
							continue;
						}
						const uint32_t ncM = M + ((o+1) << shift);
						if (lo->first >= ncM)
						{
							continue;
						}
						// use only non-empty nodes
						todo.emplace(node_info(lo, clevel));
					}
				}
			}
		}

		void get_clashing_pts_greedy(
			const Canonic* c,
			std::vector<point3r>& cpts) const
		{
			cpts.clear();
			auto iend = ids.cend();
			auto nend = bst.cend();
			{for (auto nit = bst.cbegin(); nit != bst.cend(); ++nit)
			{
				const point3i lpt = dec(nit->first);
				if (c->clashing(mapper.loc2glo(lpt), leaf_diam))
				{
					auto nnit = nit;
					auto ciend = (++nnit != nend) ? nnit->second : iend;
					{for (auto it = nit->second; it != ciend; ++it)
					{
						cpts.push_back(pts.at(*it));
					}}
				}
			}}
		}

	};


} /* namespace shade */

#endif /* ztree_h */
