/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

// #define GMP_MPN	true	// 64-bit GMP is required. Low-level functions are implemented using GMP: on Windows, mpn has limit of 2^(31 + 6) bits (41 billion digits).

#include <gmp.h>

#ifndef GMP_MPN

#ifndef __x86_64
inline uint64_t _addc(const uint64_t x, const uint64_t y, uint64_t & carry)
{
	const __uint128_t t = x + __uint128_t(y) + carry;
	carry = uint64_t(t >> 64);
	return uint64_t(t);
}

inline uint64_t _subb(const uint64_t x, const uint64_t y, int64_t & borrow)
{
	const __int128_t t = x - __int128_t(y) + borrow;
	borrow = int64_t(t >> 64);
	return uint64_t(t);
}
#endif

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
	const uint64_t x_0 = x[0] + n; x[0] = x_0;
	if (x_0 >= n) return 0;
	for (size_t i = 1; i < size; ++i)	// carry
	{
		const uint64_t x_i = x[i] + 1; x[i] = x_i;
		if (x_i != 0) return 0;
	}
	return 1;
#endif
}

// size > 0
inline void g_sub_1(uint64_t * const x, const size_t size, const uint64_t n)
{
#ifdef GMP_MPN
	mpn_sub_1(mp_ptr(x), mp_srcptr(x), mp_size_t(size), mp_limb_t(n));
#else
	const uint64_t x_0 = x[0]; x[0] = x_0 - n;
	if (x_0 >= n) return;
	for (size_t i = 1; i < size; ++i)	// borrow
	{
		const uint64_t x_i = x[i]; x[i] = x_i - 1;
		if (x_i != 0) return;
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
#ifdef __x86_64
	unsigned char carry = 0;
	asm volatile
	(
		"movq	%[y_size], %%rcx\n\t"
		"movq	%[x], %%rsi\n\t"
		"movq	%[y], %%rbx\n\t"
		"movq	%[z], %%rdi\n\t"
		"movq	%%rcx, %%rdx\n\t"
		"shrq	$2, %%rcx\n\t"
		"andq	$3, %%rdx\n\t"
		"incq	%%rdx\n\t"
		"testq	%%rcx, %%rcx\n\t"
		"clc\n\t"

		"jz		enda%=\n\t"

		"loopa%=:\n\t"
		"movq	(%%rsi), %%rax\n\t"
		"adcq	(%%rbx), %%rax\n\t"
		"movq	%%rax, (%%rdi)\n\t"
		"movq	8(%%rsi), %%rax\n\t"
		"adcq	8(%%rbx), %%rax\n\t"
		"movq	%%rax, 8(%%rdi)\n\t"
		"movq	16(%%rsi), %%rax\n\t"
		"adcq	16(%%rbx), %%rax\n\t"
		"movq	%%rax, 16(%%rdi)\n\t"
		"movq	24(%%rsi), %%rax\n\t"
		"adcq	24(%%rbx), %%rax\n\t"
		"movq	%%rax, 24(%%rdi)\n\t"

		"leaq	32(%%rsi), %%rsi\n\t"
		"leaq	32(%%rbx), %%rbx\n\t"
		"leaq	32(%%rdi), %%rdi\n\t"
		"decq	%%rcx\n\t"
		"jnz	loopa%=\n\t"
		"enda%=:\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	(%%rsi), %%rax\n\t"
		"adcq	(%%rbx), %%rax\n\t"
		"movq	%%rax, (%%rdi)\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	8(%%rsi), %%rax\n\t"
		"adcq	8(%%rbx), %%rax\n\t"
		"movq	%%rax, 8(%%rdi)\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	16(%%rsi), %%rax\n\t"
		"adcq	16(%%rbx), %%rax\n\t"
		"movq	%%rax, 16(%%rdi)\n\t"

		"endb%=:\n\t"
		"setc	%[carry]\n\t"

		: [carry] "=rm" (carry)
		: [x] "rm" (x), [y] "rm" (y), [z] "rm" (z), [y_size] "rm" (y_size)
		: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "cc", "memory"
	);
#else
	uint64_t carry = 0;
	for (size_t i = 0; i < y_size; ++i) z[i] = _addc(x[i], y[i], carry);
#endif
	size_t i = y_size;
	for (; i < x_size; ++i)
	{
		if (carry == 0) break;
		const uint64_t z_i = x[i] + 1; z[i] = z_i;
		carry = (z_i == 0) ? 1 : 0;
	}
	if (z != x) { for (; i < x_size; ++i) z[i] = x[i]; }
	return uint64_t(carry);
#endif
}

// x_size >= y_size > 0
inline void g_sub(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
#ifdef GMP_MPN
	mpn_sub(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
#else
#ifdef __x86_64
	unsigned char borrow = 0;
	asm volatile
	(
		"movq	%[y_size], %%rcx\n\t"
		"movq	%[x], %%rsi\n\t"
		"movq	%[y], %%rbx\n\t"
		"movq	%[z], %%rdi\n\t"
		"movq	%%rcx, %%rdx\n\t"
		"shrq	$2, %%rcx\n\t"
		"andq	$3, %%rdx\n\t"
		"incq	%%rdx\n\t"
		"testq	%%rcx, %%rcx\n\t"
		"clc\n\t"

		"jz		enda%=\n\t"

		"loopa%=:\n\t"
		"movq	(%%rsi), %%rax\n\t"
		"sbbq	(%%rbx), %%rax\n\t"
		"movq	%%rax, (%%rdi)\n\t"
		"movq	8(%%rsi), %%rax\n\t"
		"sbbq	8(%%rbx), %%rax\n\t"
		"movq	%%rax, 8(%%rdi)\n\t"
		"movq	16(%%rsi), %%rax\n\t"
		"sbbq	16(%%rbx), %%rax\n\t"
		"movq	%%rax, 16(%%rdi)\n\t"
		"movq	24(%%rsi), %%rax\n\t"
		"sbbq	24(%%rbx), %%rax\n\t"
		"movq	%%rax, 24(%%rdi)\n\t"

		"leaq	32(%%rsi), %%rsi\n\t"
		"leaq	32(%%rbx), %%rbx\n\t"
		"leaq	32(%%rdi), %%rdi\n\t"
		"decq	%%rcx\n\t"
		"jnz	loopa%=\n\t"
		"enda%=:\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	(%%rsi), %%rax\n\t"
		"sbbq	(%%rbx), %%rax\n\t"
		"movq	%%rax, (%%rdi)\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	8(%%rsi), %%rax\n\t"
		"sbbq	8(%%rbx), %%rax\n\t"
		"movq	%%rax, 8(%%rdi)\n\t"

		"decq	%%rdx\n\t"
		"jz		endb%=\n\t"

		"movq	16(%%rsi), %%rax\n\t"
		"sbbq	16(%%rbx), %%rax\n\t"
		"movq	%%rax, 16(%%rdi)\n\t"

		"endb%=:\n\t"
		"setc	%[borrow]\n\t"

		: [borrow] "=rm" (borrow)
		: [x] "rm" (x), [y] "rm" (y), [z] "rm" (z), [y_size] "rm" (y_size)
		: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "cc", "memory"
	);
#else
	int64_t borrow = 0;
	for (size_t i = 0;  i < y_size; ++i) z[i] = _subb(x[i], y[i], borrow);
#endif
	size_t i = y_size;
	for (; i < x_size; ++i)
	{
		if (borrow == 0) break;
		const uint64_t x_i = x[i]; z[i] = x_i - 1;
		borrow = (x_i == 0) ? 1 : 0;
	}
	if (z != x) { for (; i < x_size; ++i) z[i] = x[i]; }
#endif
}

// x_size >= y_size > 0, z_size = x_size + y_size, z != x, z != y
inline void g_mul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	mpn_mul(mp_ptr(z), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
}

// size > 0, z_size = 2 * x_size, z != x
inline void g_mul(uint64_t * const z, const uint64_t * const x, const uint64_t * const y, const size_t size)
{
	mpn_mul_n(mp_ptr(z), mp_srcptr(x), mp_srcptr(y), mp_size_t(size));
}

// size > 0, z_size = 2 * x_size, z != x
inline void g_sqr(uint64_t * const z, const uint64_t * const x, const size_t size)
{
	mpn_sqr(mp_ptr(z), mp_srcptr(x), mp_size_t(size));
}

inline void g_get_str(char * const str, const uint64_t * const x, const size_t size)
{
	const size_t str_size = mpn_get_str((unsigned char *)str, 10, mp_ptr(x), mp_size_t(size));
	for (size_t i = 0; i < str_size; ++i) str[i] += '0';
	str[str_size] = '\0';
}
