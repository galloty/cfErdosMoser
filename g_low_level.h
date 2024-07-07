/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

#include <gmp.h>

// Low-level functions implemented using GMP. But on Windows, mpn has limit of 2^(31 + 6) bits (41 billion digits).

inline void g_zero(uint64_t * const x, const size_t size)
{
	for (size_t i = 0; i < size; ++i) x[i] = 0;
}

inline void g_copy(uint64_t * const y, const uint64_t * const x, const size_t size)
{
	for (size_t i = 0; i < size; ++i) y[i] = x[i];
}

inline void g_copy_rev(uint64_t * const y, const uint64_t * const x, const size_t size)
{
	for (size_t i = 0, j = size - 1; i < size; ++i, --j) y[j] = x[j];
}

inline int g_cmp(const uint64_t * const x, const uint64_t * const y, const size_t size)
{
	for (size_t i = 0, j = size - 1; i < size; ++i, --j)
	{
		const uint64_t x_j = x[j], y_j = y[j];
		if (x_j > y_j) return 1; else if (x_j < y_j) return -1;
	}
	return 0;
}

inline uint64_t g_add_1(uint64_t * const y, const uint64_t * const x, const size_t size, const uint64_t n)
{
	return (uint64_t)mpn_add_1(mp_ptr(y), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
}

inline void g_sub_1(uint64_t * const y, const uint64_t * const x, const size_t size, const uint64_t n)
{
	mpn_sub_1(mp_ptr(y), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
}

inline uint64_t g_mul_1(uint64_t * const y, const uint64_t * const x, const size_t size, const uint64_t n)
{
	return (uint64_t)mpn_mul_1(mp_ptr(y), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
}

inline uint64_t g_mod_1(const uint64_t * const x, const size_t size, const uint64_t n)
{
	return (uint64_t)mpn_mod_1(mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
}

inline uint64_t g_add(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	// _addcarry_u64 
	return (uint64_t)mpn_add(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

inline uint64_t g_sub(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	return (uint64_t)mpn_sub(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

inline uint64_t g_mul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	return (uint64_t)mpn_mul(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

inline void g_tdiv_qr(uint64_t * const q, uint64_t * const r, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	mpn_tdiv_qr(mp_ptr(q), mp_ptr(r), mp_size_t(0), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

inline void g_get_str(char * const str, const uint64_t * const x, const size_t size)
{
	const size_t str_size = mpn_get_str((unsigned char *)str, 10, mp_ptr(x), mp_size_t(size));
	for (size_t i = 0; i < str_size; ++i) str[i] += '0';
	str[str_size] = '\0';
}
