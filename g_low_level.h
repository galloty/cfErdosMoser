/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <immintrin.h>

#ifdef GMP_MPN

// Low-level functions are implemented using GMP. On Windows, mpn has limit of 2^(31 + 6) bits (41 billion digits).
#include <gmp.h>

#else

#include "mod64.h"

inline uint64_t _addc(const uint64_t x, const uint64_t y, uint64_t & carry)
{
	const __uint128_t t = x + __uint128_t(y) + carry;
	carry = uint64_t(t >> 64);
	return uint64_t(t);
}

inline uint64_t _subb(const uint64_t x, const uint64_t y, uint64_t & borrow)
{
	const __uint128_t t = x - __uint128_t(y) - borrow;
	borrow = uint64_t(t >> 64) & 1u;
	return uint64_t(t);
}

inline uint64_t _mulc(const uint64_t x, const uint64_t y, uint64_t & carry)
{
	const __uint128_t t = x * __uint128_t(y) + carry;
	carry = uint64_t(t >> 64);
	return uint64_t(t);
}

// Modular SIMD arithmetic in Mathemagix, Joris van der Hoeven, Grégoire Lecerf, Guillaume Quintin, Chapter 2.2.

inline int _mod_exp(const uint64_t p) { int r = 0; while ((p >> r) != 0) ++r; return r; }	// 2^{r-1} <= p < 2^r

inline uint64_t _mod_invert(const uint64_t p, const int r)
{
	const int n = 62, s = r - 2, t = n + 1;
	return uint64_t((__uint128_t(1) << (s + t)) / p);
}

inline uint64_t _mod_62(const __uint128_t a, const uint64_t p, const uint64_t p_inv, const int r)
{
	const int n = 62, s = r - 2, t = n + 1;

	// a < 2^{2r}, s = r - 2 => b < 2^{r + 2} < 2^64
	const uint64_t b = uint64_t(a >> s);
	// p_inv < 2^{r - 2 + n + 1 - (r - 1)} = 2^n = 2^62 => c < 2^{64 + 62 - (62 + 1)} = 2^63
	const uint64_t c = uint64_t((b * __uint128_t(p_inv)) >> t);
	const uint64_t d = uint64_t(a - __uint128_t(c) * p);
	// The number of iteration is 1 (see Table 3)
	return (d >= p) ? d - p : d;
}

#endif

#include <gmp.h>	// TODO remove

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
#ifdef GMP_MPN
	return mpn_cmp(mp_srcptr(x), mp_srcptr(y), mp_size_t(n));
#else
	for (size_t i = 0, j = size - 1; i < size; ++i, --j)
	{
		const uint64_t x_j = x[j], y_j = y[j];
		if (x_j > y_j) return 1; else if (x_j < y_j) return -1;
	}
	return 0;
#endif
}

// size > 0
inline uint64_t g_add_1(uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_add_1(mp_ptr(x), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t carry = 0;
	x[0] = _addc(x[0], n, carry);
	for (size_t i = 1; i < size; ++i)
	{
		if (carry == 0) break;
		x[i] = _addc(x[i], 0, carry);
	}
	return carry;
#endif
}

// size > 0
inline void g_sub_1(uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	mpn_sub_1(mp_ptr(x), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t borrow = 0;
	x[0] = _subb(x[0], n, borrow);
	for (size_t i = 1; i < size; ++i)
	{
		if (borrow == 0) break;
		x[i] = _subb(x[i], 0, borrow);
	}
#endif
}

// size > 0
inline uint64_t g_mul_1(uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_mul_1(mp_ptr(x), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t carry = 0;
	for (size_t i = 0; i < size; ++i) x[i] = _mulc(x[i], n, carry);
	return carry;
#endif
}

inline void g_mod_init(const uint64_t p, uint64_t & p_inv, int & e)
{
	e = _mod_exp(p);
	p_inv = _mod_invert(p, e);
}

// size > 0, n < 2^62
inline uint64_t g_mod_1(const uint64_t * const x, const size_t size, const uint64_t n, const uint64_t n_inv, const int e, const uint64_t f)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_mod_1(mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t remainder = 0;
	for (size_t i = 0, j = size - 1; i < size; ++i, --j)
	{
		// Compute remainder * (2^64 mod n) rather than remainder * 2^64 such that t < n^2
		const __uint128_t t = remainder * __uint128_t(f) + x[j];
		remainder = _mod_62(t, n, n_inv, e);
	}

	return remainder;
#endif
}

// x_size >= y_size > 0
inline uint64_t g_add(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_add(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
#else
	uint64_t carry = 0;
	for (size_t i = 0; i < y_size; ++i) z[i] = _addc(x[i], y[i], carry);
	for (size_t i = y_size; i < x_size; ++i) z[i] = _addc(x[i], 0, carry);
	return carry;
#endif
}

inline void g_sub(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
#ifdef GMP_MPN
	mpn_sub(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
#else
	uint64_t borrow = 0;
	for (size_t i = 0; i < y_size; ++i) z[i] = _subb(x[i], y[i], borrow);
	for (size_t i = y_size; i < x_size; ++i) z[i] = _subb(x[i], 0, borrow);
#endif
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
