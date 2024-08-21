/*
Copyright 2024, Yves Gallot

SSGA is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#define _USE_MATH_DEFINES
#include <cmath>

#include <gmp.h>
#include <omp.h>

#include "heap.h"

#define finline	__attribute__((always_inline))
#ifdef __x86_64
	#define x64_asm
#endif

// A vector of l elements, operation are modulo 2^n + 1
class ModVector
{
private:
	const size_t _size, _n, _l;
	uint64_t * const _buf;
	uint64_t * const _d;

#ifndef x64_asm
	finline static uint64_t _addc(const uint64_t x, const uint64_t y, uint64_t & carry)
	{
		const __uint128_t t = x + __uint128_t(y) + carry;
		carry = uint64_t(t >> 64);
		return uint64_t(t);
	}

	finline static uint64_t _subb(const uint64_t x, const uint64_t y, int64_t & borrow)
	{
		const __int128_t t = x - __int128_t(y) + borrow;
		borrow = int64_t(t >> 64);
		return uint64_t(t);
	}
#endif

	finline static void _get(const uint64_t * const x, const size_t x_size, uint64_t * const dst, size_t & dst_size)
	{
		size_t size = x_size; while ((size != 0) && (x[size - 1] == 0)) --size;
		for (size_t i = 0; i < size; ++i) dst[i] = x[i];
		dst_size = size;
	}

	finline static void _set(uint64_t * const x, const size_t x_size, const uint64_t * const src, const size_t src_size)
	{
		for (size_t i = 0; i < src_size; ++i) x[i] = src[i];
		for (size_t i = src_size; i < x_size; ++i) x[i] = 0;
	}

	// x += a
	finline static bool _add(uint64_t * const x, const size_t size, const uint64_t a)
	{
		const uint64_t x_0 = x[0] + a; x[0] = x_0;
		if (x_0 >= a) return false;
		for (size_t i = 1; i < size; ++i)	// carry
		{
			const uint64_t x_i = x[i] + 1; x[i] = x_i;
			if (x_i != 0) return false;
		}
		return true;
	}

	// x -= a
	finline static bool _sub(uint64_t * const x, const size_t size, const uint64_t a)
	{
		const uint64_t x_0 = x[0]; x[0] = x_0 - a;
		if (x_0 >= a) return false;
		for (size_t i = 1; i < size; ++i)	// borrow
		{
			const uint64_t x_i = x[i]; x[i] = x_i - 1;
			if (x_i != 0) return false;
		}
		return true;
	}

	// x' = x + y, y' = x - y
	finline static bool _add_sub(uint64_t * const x, uint64_t * const y, const size_t size)
	{
#ifdef x64_asm
		const size_t size_4 = size / 4;		// size = 4 * size_4 + 1
		char borrow = 0;
		asm volatile
		(
			"movq	%[size_4], %%rcx\n\t"
			"movq	%[x], %%rsi\n\t"
			"movq	%[y], %%rdi\n\t"
			"xorq	%%rax, %%rax\n\t"		// carry of sbb
			"xorq	%%rdx, %%rdx\n\t"		// carry of adc
			"clc\n\t"

			"loop%=:\n\t"
			"neg	%%al\n\t"				// CF is set to 0 if dl = 0, otherwise it is set to 1

			"movq	(%%rsi), %%rbx\n\t"
			"movq	(%%rdi), %%r9\n\t"
			"movq	%%rbx, %%r8\n\t"		// r8 = x[0], r9 = y[0]
			"sbbq	%%r9, %%rbx\n\t"
			"movq	%%rbx, (%%rdi)\n\t"		// y[i] = rbx = x[i] - y[i]

			"movq	8(%%rsi), %%rbx\n\t"
			"movq	8(%%rdi), %%r11\n\t"
			"movq	%%rbx, %%r10\n\t"		// r10 = x[1], r11 = y[1]
			"sbbq	%%r11, %%rbx\n\t"
			"movq	%%rbx, 8(%%rdi)\n\t"

			"movq	16(%%rsi), %%rbx\n\t"
			"movq	16(%%rdi), %%r13\n\t"
			"movq	%%rbx, %%r12\n\t"		// r12 = x[2], r13 = y[2]
			"sbbq	%%r13, %%rbx\n\t"
			"movq	%%rbx, 16(%%rdi)\n\t"

			"movq	24(%%rsi), %%rbx\n\t"
			"movq	24(%%rdi), %%r15\n\t"
			"movq	%%rbx, %%r14\n\t"		// r14 = x[3], r15 = y[3]
			"sbbq	%%r15, %%rbx\n\t"
			"movq	%%rbx, 24(%%rdi)\n\t"

			"setc	%%al\n\t"
			"neg	%%dl\n\t"

			"adcq	%%r9, %%r8\n\t"
			"movq	%%r8, (%%rsi)\n\t"		// x[i] = x[i] + y[i]
			"adcq	%%r11, %%r10\n\t"
			"movq	%%r10, 8(%%rsi)\n\t"
			"adcq	%%r13, %%r12\n\t"
			"movq	%%r12, 16(%%rsi)\n\t"
			"adcq	%%r15, %%r14\n\t"
			"movq	%%r14, 24(%%rsi)\n\t"

			"setc	%%dl\n\t"

			"decq	%%rcx\n\t"
			"addq	$32, %%rsi\n\t"
			"addq	$32, %%rdi\n\t"
			"testq	%%rcx, %%rcx\n\t"
			"jne	loop%=\n\t"

			"neg	%%al\n\t"

			"movq	(%%rsi), %%rbx\n\t"
			"movq	(%%rdi), %%r9\n\t"
			"movq	%%rbx, %%r8\n\t"
			"sbbq	%%r9, %%rbx\n\t"
			"movq	%%rbx, (%%rdi)\n\t"

			"setc	%[borrow]\n\t"
			"neg	%%dl\n\t"

			"adcq	%%r9, %%r8\n\t"
			"movq	%%r8, (%%rsi)\n\t"

			: [borrow] "=rm" (borrow)
			: [x] "rm" (x), [y] "rm" (y), [size_4] "rm" (size_4)
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory"
		);
#else
	uint64_t carry = 0; int64_t borrow = 0;
	for (size_t i = 0; i < size; ++i)
	{
		const uint64_t x_i = x[i], y_i = y[i];
		x[i] = _addc(x_i, y_i, carry);
		y[i] = _subb(x_i, y_i, borrow);
	}
#endif
		return (borrow != 0);
	}

	// y = x << s, 0 <= s < 64
	finline static void _lshift(uint64_t * const y, const uint64_t * const x, const size_t size, const unsigned int s)
	{
		if (s > 0)
		{
			uint64_t prev = x[0]; y[0] = prev << s;
			for (size_t i = 1; i < size; ++i)
			{
				const uint64_t x_i = x[i]; y[i] = (x_i << s) | (prev >> (64 - s)); prev = x_i;
			}
		}
		else for (size_t i = 0; i < size; ++i) y[i] = x[i];
	}

	// x ?>= 2^n + 1
	finline static bool _ge_F(const uint64_t * const x, const size_t size)
	{
		if (x[size - 1] != 1) return (x[size - 1] > 1);
		for (size_t i = 0; i < size - 1; ++i) if (x[i] != 0) return true;
		return false;
	}

	// x += 2^n + 1
	finline static void _add_F(uint64_t * const x, const size_t size) { x[size - 1] += 1; _add(x, size, 1); }

	// x -= 2^n + 1
	finline static void _sub_F(uint64_t * const x, const size_t size) { x[size - 1] -= 1; _sub(x, size, 1); }

	// x = y (mod 2^n + 1)
	finline static void _mod_F(uint64_t * const x, const size_t size, const uint64_t * const y)
	{
#ifdef x64_asm
		char borrow = 0;
		asm volatile
		(
			"movq	%[size], %%rcx\n\t"
			"movq	%[y], %%rsi\n\t"
			"movq	%[x], %%rdi\n\t"
			"leaq	-8(%%rsi,%%rcx,8), %%rbx\n\t"	// y[size - 1]
			"shrq	$2,	%%rcx\n\t"					// size = 4 * size_4 + 1
			"clc\n\t"

			"loop%=:\n\t"
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
			"jnz	loop%=\n\t"

			"movq	$0, %%rax\n\t"
			"movq	(%%rbx), %%rdx\n\t"
			"sbbq	%%rdx, %%rax\n\t"
			"movq	%%rax, (%%rdi)\n\t"

			"setc	%[borrow]\n\t"

			: [borrow] "=rm" (borrow)
			: [x] "rm" (x), [y] "rm" (y), [size] "rm" (size)
			: "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "cc", "memory"
		);
#else
		// y mod 2^n - y / 2^n
	 	int64_t borrow = 0;
	 	for (size_t i = 0; i < size - 1; ++i) x[i] = _subb(y[i], y[size - 1 + i], borrow);
		x[size - 1] = _subb(0, y[2 * size - 2], borrow);
#endif
	 	if (borrow != 0) _add_F(x, size);
	}

	// y = x << (64 * s_64) (mod 2^n + 1), s_64 < x_size - 1
	finline static void _lshift_mod_F(uint64_t * const y, const uint64_t * const x, const size_t size, const size_t s_64)
	{
		for (size_t i = 0; i < s_64; ++i) y[i] = ~x[size - 1 - s_64 + i];
		for (size_t i = s_64; i < size - 1; ++i) y[i] = x[i - s_64];
		y[size - 1] = 0;

		if (s_64 != 0)
		{
			const bool carry = _add(y, s_64, 1);	// two's complement
			// If carry then we have 0 = -0 otherwise sub 1 to next term
			if (!carry) { const bool borrow = _sub(&y[s_64], size - s_64, 1); if (borrow) _add_F(y, size); }
		}

		const bool borrow = _sub(&y[s_64], size - s_64, x[size - 1]); if (borrow) _add_F(y, size);
	}

	// y = -x << (64 * s_64) (mod 2^n + 1), s_64 < x_size - 1
	finline static void _nlshift_mod_F(uint64_t * const y, const uint64_t * const x, const size_t size, const size_t s_64)
	{
	 	for (size_t i = 0; i < s_64; ++i) y[i] = x[size - 1 - s_64 + i];
	 	for (size_t i = s_64; i < size - 1; ++i) y[i] = ~x[i - s_64];
		y[size - 1] = 0;

		const bool carry = _add(&y[s_64], size - 1 - s_64, 1);	// two's complement
		// If carry then we have 0 = -0 otherwise add 1 = -2^n
		if (!carry) _add(y, size, 1);

		_add(&y[s_64], size - s_64, x[size - 1]);
		if (_ge_F(y, size)) _sub_F(y, size);
	}

	void add_sub(const size_t i, const size_t j)
	{
		const size_t size = _size;
		uint64_t * const d_i = &_d[i * (size + _gap)];
		uint64_t * const d_j = &_d[j * (size + _gap)];
		const bool borrow = _add_sub(d_i, d_j, size);
		if (_ge_F(d_i, size)) _sub_F(d_i, size);
		if (borrow) _add_F(d_j, size);
	}

	void lshift(const size_t i, const size_t s, uint64_t * const buf)
	{
		const size_t size = _size;
		uint64_t * const d_i = &_d[i * (size + _gap)];
		_lshift(buf, d_i, size, s % 64);
		_lshift_mod_F(d_i, buf, size, s / 64);
	}

	void rshift(const size_t i, const size_t s, uint64_t * const buf)
	{
		// 2^n = -1 then 2^-s = -2^{n - s}
		const size_t size = _size, ls = _n - s;
		uint64_t * const d_i = &_d[i * (size + _gap)];
		_lshift(buf, d_i, size, ls % 64);
		_nlshift_mod_F(d_i, buf, size, ls / 64);
	}

	void negacyclic(const ModVector & rhs, const size_t i, uint64_t * const buf)
	{
		const size_t size = _size;
		uint64_t * const d_i = &_d[i * (size + _gap)];
		mpn_mul_n(mp_ptr(buf), mp_srcptr(d_i), mp_srcptr(&rhs._d[i * (size + _gap)]), mp_size_t(size));
		_mod_F(d_i, size, buf);
	}

	static const size_t _gap = 7;	// Cache line size is 64 bytes

public:
	ModVector(const size_t n, const size_t l) : _size(n / 64 + 1), _n(n), _l(l),
		_buf(static_cast<uint64_t *>(Heap::get_instance().aligned_alloc(4 * 2 * _size * sizeof(uint64_t)))),
		_d(static_cast<uint64_t *>(Heap::get_instance().aligned_alloc(l * (_size + _gap) * sizeof(uint64_t)))) {}
	virtual ~ModVector()
	{
		Heap & heap = Heap::get_instance();
		heap.aligned_free(_buf, 4 * 2 * _size * sizeof(uint64_t));
		heap.aligned_free(_d, _l * (_size + _gap) * sizeof(uint64_t));
	}

	size_t get_size() const { return _size; }

	void get(const size_t i, uint64_t * const x, size_t & x_size) const { _get(&_d[i * (_size + _gap)], _size, x, x_size); }
	void set(const size_t i, const uint64_t * const x, const size_t x_size) { _set(&_d[i * (_size + _gap)], _size, x, x_size); }

	void mul(const ModVector & rhs, const size_t m, const size_t j, const size_t e, uint64_t * const buf)
	{
		// previous root is r^2 = 2^e, new root is r = 2^{e/2}
		const size_t e_2 = e / 2;

		for (size_t i = 0; i < m; ++i) { lshift(j + i + 1 * m, e_2, buf); add_sub(j + i + 0 * m, j + i + 1 * m); }

		if (m > 1)
		{
			mul(rhs, m / 2, j + 0 * m, e_2, buf);		//  r = 2^{e/2}
			mul(rhs, m / 2, j + 1 * m, e_2 + _n, buf);	// -r = 2^{e/2 + n}
		}
		else
		{
			negacyclic(rhs, j + 0 * m, buf);
			negacyclic(rhs, j + 1 * m, buf);
		}

		for (size_t i = 0; i < m; ++i) { add_sub(j + i + 0 * m, j + i + 1 * m); rshift(j + i + 1 * m, e_2, buf); }
	}

	void mul_Mersenne(const ModVector & rhs, const size_t m, const size_t j)
	{
		// We have e = 0, root is 1
		for (size_t i = 0; i < m; ++i) add_sub(j + i + 0 * m, j + i + 1 * m);

		if (m > 1)
		{
			mul_Mersenne(rhs, m / 2, j + 0 * m);	// root of 1 is still 1
			mul(rhs, m / 2, j + 1 * m, _n, _buf);	// root is -1 = 2^n
		}
		else
		{
			negacyclic(rhs, j + 0 * m, _buf);
			negacyclic(rhs, j + 1 * m, _buf);
		}

		for (size_t i = 0; i < m; ++i) add_sub(j + i + 0 * m, j + i + 1 * m);
	}

	void mul_Mersenne_0(const ModVector & rhs, const size_t m, const size_t k, const int nthreads)
	{
		const size_t m_2 = m / 2, n_2 = _n / 2;

		for (size_t i = 0; i < m_2; ++i)
		{
			add_sub(i + 0 * m_2, i + 2 * m_2); add_sub(i + 1 * m_2, i + 3 * m_2);
			add_sub(i + 0 * m_2, i + 1 * m_2); lshift(i + 3 * m_2, n_2, _buf); add_sub(i + 2 * m_2, i + 3 * m_2);
		}

		if (nthreads == 4)
		{
#pragma omp parallel
		{
			const unsigned int thread_id = static_cast<unsigned int>(omp_get_thread_num());
			if      (thread_id == 0) mul_Mersenne(rhs, m_2 / 2, 0 * m_2);
			else if (thread_id == 1) mul(rhs, m_2 / 2, 1 * m_2, 2 * n_2, &_buf[1 * 2 * _size]);
			else if (thread_id == 2) mul(rhs, m_2 / 2, 2 * m_2, 1 * n_2, &_buf[2 * 2 * _size]);
			else if (thread_id == 3) mul(rhs, m_2 / 2, 3 * m_2, 3 * n_2, &_buf[3 * 2 * _size]);
		}
		}
		else
		{
			mul_Mersenne(rhs, m_2 / 2, 0 * m_2);
			mul(rhs, m_2 / 2, 1 * m_2, 2 * n_2, _buf);
			mul(rhs, m_2 / 2, 2 * m_2, 1 * n_2, _buf);
			mul(rhs, m_2 / 2, 3 * m_2, 3 * n_2, _buf);

		}

		// Components are not halved during the reverse transform then multiply outputs by 1 / 2^k
		for (size_t i = 0; i < m_2; ++i)
		{
			add_sub(i + 0 * m_2, i + 1 * m_2); add_sub(i + 2 * m_2, i + 3 * m_2); rshift(i + 3 * m_2, n_2, _buf);
			add_sub(i + 0 * m_2, i + 2 * m_2); rshift(i + 0 * m_2, k, _buf); rshift(i + 2 * m_2, k, _buf);
			add_sub(i + 1 * m_2, i + 3 * m_2); rshift(i + 1 * m_2, k, _buf); rshift(i + 3 * m_2, k, _buf);
		} 
	}

	void forward(const size_t m, const size_t j, const size_t e, uint64_t * const buf)
	{
		const size_t e_2 = e / 2;
		for (size_t i = 0; i < m; ++i) { lshift(j + i + 1 * m, e_2, buf); add_sub(j + i + 0 * m, j + i + 1 * m); }
		if (m > 1)
		{
			forward(m / 2, j + 0 * m, e_2, buf);
			forward(m / 2, j + 1 * m, e_2 + _n, buf);
		}
	}

	void forward_Mersenne(const size_t m, const size_t j)
	{
		for (size_t i = 0; i < m; ++i) add_sub(j + i + 0 * m, j + i + 1 * m);
		if (m > 1)
		{
			forward_Mersenne(m / 2, j + 0 * m);
			forward(m / 2, j + 1 * m, _n, _buf);
		}
	}

	void forward_Mersenne_0(const size_t m, const int nthreads)
	{
		const size_t m_2 = m / 2, n_2 = _n / 2;

		for (size_t i = 0; i < m_2; ++i)
		{
			add_sub(i + 0 * m_2, i + 2 * m_2); add_sub(i + 1 * m_2, i + 3 * m_2);
			add_sub(i + 0 * m_2, i + 1 * m_2); lshift(i + 3 * m_2, n_2, _buf); add_sub(i + 2 * m_2, i + 3 * m_2);
		}

		if (nthreads == 4)
		{
#pragma omp parallel
		{
			const unsigned int thread_id = static_cast<unsigned int>(omp_get_thread_num());
			if      (thread_id == 0) forward_Mersenne(m_2 / 2, 0 * m_2);
			else if (thread_id == 1) forward(m_2 / 2, 1 * m_2, 2 * n_2, &_buf[1 * 2 * _size]);
			else if (thread_id == 2) forward(m_2 / 2, 2 * m_2, 1 * n_2, &_buf[2 * 2 * _size]);
			else if (thread_id == 3) forward(m_2 / 2, 3 * m_2, 3 * n_2, &_buf[3 * 2 * _size]);
		}
		}
		else
		{
			forward_Mersenne(m_2 / 2, 0 * m_2);
			forward(m_2 / 2, 1 * m_2, 2 * n_2, _buf);
			forward(m_2 / 2, 2 * m_2, 1 * n_2, _buf);
			forward(m_2 / 2, 3 * m_2, 3 * n_2, _buf);
		}
	}
};

class SSG
{
private:
	const unsigned int _k;
	const size_t _M, _n, _l;
	ModVector _x, _y;
	static int _nthreads;

	// v is a vector of l M-bit slices of x
	void set_vector(ModVector & v, const uint64_t * const x, const size_t size)
	{
		const size_t M_64 = _M / 64, M_mask = (size_t(1) << (_M % 64)) - 1;

		uint64_t * const r = new uint64_t[M_64 + 2];

		for (size_t j = 0, l = _l; j < l; ++j)
		{
			const size_t bit_index = j * _M, index = bit_index / 64;

			if (index < size)
			{
				const size_t left = size - index;
				for (size_t i = 0, n = std::min(left, M_64 + 2); i < n; ++i) r[i] = x[index + i];
				for (size_t i = left; i < M_64 + 2; ++i) r[i] = 0;

				const unsigned int s = bit_index % 64;
				if (s != 0)	// right shift
				{
					uint64_t next = r[1]; r[0] = (r[0] >> s) | (next << (64 - s));
					for (size_t i = 1; i < M_64 + 1; ++i)
					{
						const uint64_t r_i = next; next = r[i + 1];
						r[i] = (r_i >> s) | (next << (64 - s));
					}
				}
				r[M_64] &= M_mask;

				v.set(j, r, M_64 + 1);
			}
			else v.set(j, r, 0);
		}

		delete[] r;
	}

	// Compute sum of the l slices
	void get_vector(uint64_t * const x, const size_t size, const ModVector & v)
	{
		size_t clear_index = 0;

		uint64_t * const r = new uint64_t[v.get_size() + 1];

		for (size_t j = 0, l = _l; j < l; ++j)
		{
			const size_t bit_index = j * _M, index = bit_index / 64; const unsigned int s = bit_index % 64;

			size_t r_size; v.get(j, r, r_size);
			if (r_size != 0)
			{
				if (s != 0)	// left shift
				{
					uint64_t prev = r[0]; r[0] = prev << s;
					for (size_t i = 1; i < r_size; ++i)
					{
						const uint64_t r_i = r[i]; r[i] = (r_i << s) | (prev >> (64 - s)); prev = r_i;
					}
					r[r_size] = prev >> (64 - s);
				}
				else r[r_size] = 0;

				const size_t prev_clear_index = clear_index; clear_index = std::min(index + r_size + 1, size);
				for (size_t i = prev_clear_index; i < clear_index; ++i) x[i] = 0;
				mpn_add_n(mp_ptr(&x[index]), mp_srcptr(&x[index]), mp_srcptr(r), mp_size_t(clear_index - index));
			}
		}

		delete[] r;
		for (size_t i = clear_index; i < size; ++i) x[i] = 0;
	}

	static double get_param(const size_t N, const unsigned int k, size_t & M, size_t & n)
	{
		// See Pierrick Gaudry, Alexander Kruppa, Paul Zimmermann.
		// A GMP-based implementation of Schönhage-Strassen's large integer multiplication algorithm.
		// ISSAC 2007, Jul 2007, Waterloo, Ontario, Canada. pp.167-174, ⟨10.1145/1277548.1277572⟩. ⟨inria-00126462v2⟩
		const size_t K = size_t(1) << k;
		M = (N % K == 0) ? N / K : N / K + 1;
		const size_t t = 2 * M + k;
		n = t; if (n % K != 0) n = (n / K + 1) * K;
		while (n % (64 * 4) != 0) n += K;
		return double(t) / n;	// efficiency
	}

public:
	SSG(const unsigned int k, const size_t M, const size_t n) : _k(k), _M(M), _n(n), _l(size_t(1) << k), _x(n, _l), _y(n, _l) {}
	virtual ~SSG() {}

	static void set_nthreads(const int nthreads)
	{
		if (nthreads > 1) omp_set_num_threads(4);
		_nthreads = 1;
		if (nthreads != 1)
		{
#pragma omp parallel
			{
#pragma omp single
				_nthreads = omp_get_num_threads();
			}
		}
	}

	void set_y(const uint64_t * const y, const size_t size) { set_vector(_y, y, size); _y.forward_Mersenne_0(_l / 2, _nthreads); }

	void get_x(uint64_t * const x, const size_t size) { get_vector(x, size, _x); }

	void mul_xy(const uint64_t * const x, const size_t size)
	{
		set_vector(_x, x, size);
		_x.mul_Mersenne_0(_y, _l / 2, _k, _nthreads);
	}

	void sqr(const uint64_t * const x, const size_t size)
	{
		set_vector(_x, x, size);
		_x.mul_Mersenne_0(_x, _l / 2, _k, _nthreads);
	}

	static void get_best_param(const size_t N, unsigned int & k, size_t & M, size_t & n)
	{
		k = 5;
		for (unsigned int i = 6; true; ++i)
		{
			size_t M_i, n_i; const double efficiency = get_param(N, i, M_i, n_i);
			const size_t K_i = size_t(1) << i;
			if (K_i > 2 * std::sqrt(M_i * K_i)) break;
			if (efficiency > 0.95) k = i;
		}
		get_param(N, k, M, n);
	}
};

int SSG::_nthreads = 4;
