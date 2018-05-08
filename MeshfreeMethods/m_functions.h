#/*! \file n_functions.h */

//@HEADER
// ************************************************************************
//
//                         M4_Technologies
//                 Copyright (2017) by Wentao Xu
// 
// M4 stands for Multiscale Mesh-based Meshfree Method
// Any questions, please contact:
// Wentao Xu   wx2151@columbia.edu
//
// ************************************************************************
//@HEADER

#ifndef _M4_FUNCTIONS_
#define _M4_FUNCTIONS_

#include "m_types.h"

namespace msl {
	
	template<class value_type>
	value_type mean(value_type* arr, int start, int end) {
		value_type sum = 0;
		for (int i = start; i <= end; i++)
			sum += arr[i];
		return sum;
	}

	template<class value_type>
	double avg(value_type* arr, int start, int end) {
		double sz = end - start + 1;
		return (double)mean(arr, start, end) / sz;
	}

	inline Vec3i equals(const Vec3i& a, const Vec3i& b) {
		Vec3i temp = vec3zero;
		if (a[0] == b[0]) temp[0] = 1;
		if (a[1] == b[1]) temp[1] = 1;
		if (a[2] == b[2]) temp[2] = 1;
		return temp;
	}

	inline uint64_t split(unsigned a) {
		uint64_t x = a & 0x1fffff;                  // Only look at the first 21 bits
		x = (x | x << 32) & 0x1f00000000ffff;       // shift left 32 bits, OR with self, AND 
		x = (x | x << 16) & 0x1f0000ff0000ff;       // shift left 32 bits, OR with self, AND 
		x = (x | x << 8) & 0x100f00f00f00f00f;      // shift left 32 bits, OR with self, AND 
		x = (x | x << 4) & 0x10c30c30c30c30c3;      // shift left 32 bits, OR with self, AND 
		x = (x | x << 2) & 0x1249249249249249;
		return x;
	}

	inline uint64_t mortonEncode(Vec3i cell) {
		uint64_t ans = 0;
		ans |= split(cell[0]) | (split(cell[1]) << 1) | (split(cell[2]) << 2);
		return ans;
	}
}
#endif // !_M4_FUNCTIONS_