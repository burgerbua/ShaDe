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
	
	// integer point
	using point3i = std::array<uint32_t, 3>;

	// double point
	using point3r = std::array<double, 3>;
}

bool operator==(
	const shade::point3i& l,
	const shade::point3i& r)
{
	return l[0] == r[0] && l[1] == r[1] && l[2] == r[2];
}

template <typename point3>
point3 operator-(
	const point3& l,
	const point3& r)
{
	return{ { l[0] - r[0], l[1] - r[1], l[2] - r[2] } };
}
template <typename point3>
point3 operator+(
	const point3& l,
	const point3& r)
{
	return{ { l[0] + r[0], l[1] + r[1], l[2] + r[2] } };
}
template <typename point3>
point3 operator*(
	const point3& l,
	const double val)
{
	return{ { l[0] * val, l[1] * val, l[2] * val } };
}
template <typename point3>
double dot(
	const point3& l,
	const point3& r)
{
	return l[0] * r[0] + l[1] * r[1] + l[2] * r[2];
}
template <typename point3>
point3 cross(
	const point3& l,
	const point3& r)
{
	return{ { l[1] * r[2] - l[2] * r[1], l[2] * r[0] - l[0] * r[2], l[0] * r[1] - l[1] * r[0] } };
}
template <typename point3>
double norm(
	const point3& v)
{
	return sqrt(dot(v, v));
}
template <typename point3>
point3 normalize(
	const point3& vec)
{
	const double len = norm(vec);
	return{ { vec[0] / len, vec[1] / len, vec[2] / len } };
}

#endif /* point_h */
