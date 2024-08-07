/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <immintrin.h>

#include "fastmul.h"

// #define GMP_MPN	true	// Must be 64-bit GMP

#include <gmp.h>

#ifdef GMP_MPN
// Low-level functions are implemented using GMP. On Windows, mpn has limit of 2^(31 + 6) bits (41 billion digits).
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
	borrow = uint64_t(t >> 64) & 1;
	return uint64_t(t);
}

inline uint64_t _mulc(const uint64_t x, const uint64_t y, uint64_t & carry)
{
	const __uint128_t t = x * __uint128_t(y) + carry;
	carry = uint64_t(t >> 64);
	return uint64_t(t);
}

inline uint64_t _madc(const uint64_t x, const uint64_t y, const uint64_t z, uint64_t & carry)
{
	const __uint128_t t = z + x * __uint128_t(y) + carry;
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
	return mpn_cmp(mp_srcptr(x), mp_srcptr(y), mp_size_t(size));
#else
	for (size_t i = 0, j = size - 1; i < size; ++i, --j)
	{
		const uint64_t x_j = x[j], y_j = y[j];
		if (x_j != y_j) return (x_j > y_j) ? 1 : -1;
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
inline uint64_t g_mul_1(uint64_t * const y, const uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_mul_1(mp_ptr(y), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t carry = 0;
	for (size_t i = 0; i < size; ++i) y[i] = _mulc(x[i], n, carry);
	return carry;
#endif
}

inline void g_mod_init(const uint64_t p, uint64_t & p_inv, int & e)
{
#ifdef GMP_MPN
	(void)p; (void)p_inv; (void)e;	// remove compiler warning: unused parameter
#else
	e = _mod_exp(p);
	p_inv = _mod_invert(p, e);
#endif
}

// size > 0, n < 2^62
inline uint64_t g_mod_1(const uint64_t * const x, const size_t size, const uint64_t n, const uint64_t n_inv, const int e, const uint64_t f)
{
#ifdef GMP_MPN
	(void)n_inv; (void)e; (void)f;	// remove compiler warning: unused parameter
	return (uint64_t)mpn_mod_1(mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t r = _mod_62(x[size - 1], n, n_inv, e);
	for (size_t i = 1, j = size - 2; i < size; ++i, --j)
	{
		// Compute r * (2^64 mod n) rather than r * 2^64 such that t < n^2
		const __uint128_t t = r * __uint128_t(f) + x[j];
		r = _mod_62(t, n, n_inv, e);
	}
	return r;
#endif
}

// size > 0
inline uint64_t g_mod_1(const uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_mod_1(mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	uint64_t r = uint64_t(x[size - 1] % n);
	for (size_t i = 1, j = size - 2; i < size; ++i, --j)
	{
		const __uint128_t t = (__uint128_t(r) << 64) | x[j];
		r = uint64_t(t % n);
	}
	return r;
#endif
}

// size > 0
inline uint64_t g_div_rem_1(uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	return (uint64_t)mpn_divrem_1(mp_ptr(x), 0, mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	const uint64_t t = x[size - 1];
	uint64_t r = uint64_t(t % n); x[size - 1] = t / n;
	for (size_t i = 1, j = size - 2; i < size; ++i, --j)
	{
		const __uint128_t t = (__uint128_t(r) << 64) | x[j];
		r = uint64_t(t % n); x[j] = uint64_t(t / n);
	}
	return r;
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

// x_size >= y_size > 0
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

// x_size >= y_size > 0, z_size = x_size + y_size
inline void g_mul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	mpn_mul(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

// size > 0, z_size = 2 * x_size
inline void g_sqr(uint64_t * const z, const uint64_t * const x, const size_t size)
{
	mpn_sqr(mp_ptr(z), mp_srcptr(x), mp_size_t(size));
}

inline void g_get_str(char * const str, const uint64_t * const x, const size_t size)
{
#ifdef GMP_MPN
	const size_t str_size = mpn_get_str((unsigned char *)str, 10, mp_ptr(x), mp_size_t(size));
	for (size_t i = 0; i < str_size; ++i) str[i] += '0';
	str[str_size] = '\0';
#else
	if (size == 0) { str[0] = '0'; str[1] = '\0'; return; }

	size_t qsize = size;
	uint64_t * const q = new uint64_t[size];
	g_copy(q, x, size);

	size_t k = 0;
	while (qsize != 0)
	{
		uint64_t r = 0;
		for (size_t i = 0, j = qsize - 1; i < qsize; ++i, --j)
		{
			const __uint128_t t = (__uint128_t(r) << 64) | q[j];
			q[j] = uint64_t(t / 10); r = uint64_t(t % 10);
		}
		str[k] = '0' + char(r); ++k;
		if (q[qsize - 1] == 0) --qsize;
	}

	str[k] = '\0';
	for (size_t i = 0; i < k / 2; ++i) std::swap(str[i], str[k - i - 1]);

	delete[] q;
#endif
}
