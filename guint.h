/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

#include "g_low_level.h"
#include "gfloat.h"

// Memory allocation
class Heap
{
private:
	static size_t _size, _max_size, _max_block_size;
	static size_t _size_gmp, _max_size_gmp, _max_block_size_gmp;

public:
	static void * allocate_function(const size_t size)
	{
		_size += size;
		_max_size = std::max(_max_size, _size);
		_max_block_size = std::max(_max_block_size, size);
		return malloc(size);
	}

	static void * reallocate_function(void * const ptr, const size_t old_size, const size_t new_size)
	{
		_size += new_size - old_size;
		_max_size = std::max(_max_size, _size);
		_max_block_size = std::max(_max_block_size, new_size);
		return realloc(ptr, new_size);
	}

	static void free_function(void * const ptr, const size_t size) { _size -= size; if (size > 0) free(ptr); }

	static void * allocate_function_gmp(size_t size)
	{
		_size_gmp += size;
		_max_size_gmp = std::max(_max_size_gmp, _size_gmp);
		_max_block_size_gmp = std::max(_max_block_size_gmp, size);
		return malloc(size);
	}

	static void * reallocate_function_gmp(void *ptr, size_t old_size, size_t new_size)
	{
		_size_gmp += new_size - old_size;
		_max_size_gmp = std::max(_max_size_gmp, _size_gmp);
		_max_block_size_gmp = std::max(_max_block_size_gmp, new_size);
		return realloc(ptr, new_size);
	}

	static void free_function_gmp(void *ptr, size_t size) { _size_gmp -= size; if (size > 0) free(ptr); }

public:
	Heap() { mp_set_memory_functions(allocate_function_gmp, reallocate_function_gmp, free_function_gmp); }
	virtual ~Heap() { mp_set_memory_functions(nullptr, nullptr, nullptr); }

	size_t get_size() const { return _size; }
	size_t get_max_size() const { return _max_size; }
	size_t get_max_block_size() const { return _max_block_size; }
	size_t get_size_gmp() const { return _size_gmp; }
	size_t get_max_size_gmp() const { return _max_size_gmp; }
	size_t get_max_block_size_gmp() const { return _max_block_size_gmp; }

	void reset_max_size() { _max_size = 0; _max_size_gmp = 0; _max_block_size = 0; _max_block_size_gmp = 0; }
};

size_t Heap::_size = 0, Heap::_max_size = 0, Heap::_max_block_size = 0;
size_t Heap::_size_gmp = 0, Heap::_max_size_gmp = 0, Heap::_max_block_size_gmp = 0;

// giant unsigned integer
class guint
{
private:
	size_t _alloc;
	size_t _size;
	uint64_t * _d;

private:
	void _allocate(const size_t size)
	{
		_alloc = (size / 1024 + 1) * 1024;
		_size = size;
		_d = (uint64_t *)Heap::allocate_function(_alloc * sizeof(uint64_t));
	}

	void _reallocate(const size_t size)
	{
		if ((size > _alloc) || (_alloc - size > 1024 * 1024))
		{
			const size_t alloc = (size / 1024 + 1) * 1024;
			_d = (uint64_t *)Heap::reallocate_function(_d, _alloc * sizeof(uint64_t), alloc * sizeof(uint64_t));
			_alloc = alloc;
		}
		_size = size;
	}

	void _free()
	{
		Heap::free_function(_d, _alloc * sizeof(uint64_t));
	}

	void _norm()
	{
		const uint64_t * const d = _d;
		size_t size = _size; while ((size > 0) && (d[size - 1] == 0)) --size;
		_size = size;
	}

	// Barrett reduction: y_size <= x_size <= 2 * y_size
	// x: dividend, y: divisor, *this: quotient, r: remainder
	// Sources and destinations must be different
	void _div(guint & r, const guint & x, const guint & y, const guint & y_inv)
	{
		const size_t xsize = x._size, ysize = y._size;

		// r = x >> (ysize - 1)
		if (xsize <= ysize - 1)
		{
			*this = 0u;
			r = x;
		}
		else
		{
			r._reallocate(xsize - (ysize - 1));
			g_copy(r._d, &x._d[ysize - 1], xsize - (ysize - 1));

			// *this = (r * y_inv) >> (ysize + 1)
			mul(r, y_inv);
			if (_size <= ysize + 1)
			{
				*this = 0u;
				r = x;
			}
			else
			{
				g_copy(_d, &_d[ysize + 1], _size - (ysize + 1));
				_reallocate(_size - (ysize + 1));

				// r = x - *this * y
				r.mul(*this, y);
				if (r.sub(x)) std::cout << "Warning: _div" << std::endl;	// if x > r then sub returns x - r

				size_t h = 0;
				while (r.cmp(y) >= 0) { r -= y; *this += 1u; ++h; if (h > 1) std::cout << "Warning: _div" << std::endl; }
			}
		}
	}

public:
	guint() { _allocate(0); }
	guint(const uint64_t n) { _allocate((n == 0) ? 0 : 1); _d[0] = n; }
	guint(const guint & rhs) { _allocate(rhs._size); g_copy(_d, rhs._d, rhs._size); }
	virtual ~guint() { _free(); }

	size_t get_word_count() const { return _size; }
	size_t get_byte_count() const { return get_word_count() * sizeof(uint64_t); }

	guint & operator=(const uint64_t n) { _reallocate((n == 0) ? 0 : 1); _d[0] = n; return *this; }
	guint & operator=(const guint & rhs) { if (&rhs != this) { _reallocate(rhs._size); g_copy(_d, rhs._d, rhs._size); } return *this; }

	guint & swap(guint & rhs) { std::swap(_alloc, rhs._alloc); std::swap(_size, rhs._size); std::swap(_d, rhs._d); return *this; }

	bool is_zero() const { return (_size == 0); }

	int cmp(const guint & rhs) const
	{
		const size_t size = _size, rsize = rhs._size;
		if (size != rsize) return (size > rsize) ? 1 : -1;
		return g_cmp(_d, rhs._d, size);
	}

	guint & operator+=(const uint64_t n)
	{
		const size_t size = _size;
		if (size == 0) *this = n;
		else
		{
			const uint64_t carry = g_add_1(_d, size, n);
			if (carry != 0)
			{
				_reallocate(size + 1);
				_d[size] = carry;
			}
		}
		return *this;
	}

	// return false if negative
	bool sub(const uint64_t n)
	{
		if (n == 0) return true;

		const size_t size = _size;
		uint64_t * const d = _d;

		if (size == 0) { d[0] = n; _size = 1; return false; }
		if ((size == 1) && (d[0] < n)) { d[0] = n - d[0]; return false; }

		g_sub_1(d, size, n);
		if (d[size - 1] == 0) --_size;
		return true;
	}

	// *this >= n is required
	guint & operator-=(const uint64_t n) { sub(n); return *this; }

	guint & operator*=(const uint64_t n)
	{
		const size_t size = _size;
		if ((size == 0) || (n == 0)) { *this = 0u; }
		else
		{
			const uint64_t carry = g_mul_1(_d, size, n);
			if (carry != 0)
			{
				_reallocate(size + 1);
				_d[size] = carry;
			}
		}
		return *this;
	}

	static void mod_init(const uint64_t n, uint64_t & n_inv, int & n_e) { g_mod_init(n, n_inv, n_e); }

	uint64_t mod(const uint64_t n, const uint64_t n_inv, const int n_e, const uint64_t n_f) const
	{
		if (n == 0) throw std::runtime_error("divide by zero");
		if (_size == 0) return 0;
		return g_mod_1(_d, _size, n, n_inv, n_e, n_f);
	}

	guint & operator+=(const guint & rhs)
	{
		const size_t size = _size;
		if (size == 0) *this = rhs;
		else
		{
			const size_t rsize = rhs._size;
			if (rsize != 0)
			{
				const size_t new_size = std::max(size, rsize) + 1;
				_reallocate(new_size);
				uint64_t * const d = _d;
				const uint64_t carry = (size >= rsize) ? g_add(d, d, size, rhs._d, rsize) : g_add(d, rhs._d, rsize, d, size);
				if (carry != 0) d[new_size - 1] = carry; else --_size;
			}
		}
		return *this;
	}

	// return false if negative
	bool sub(const guint & rhs)
	{
		const size_t size = _size, rsize = rhs._size;

		if (rsize == 0) return true;
		if (size == 0) { *this = rhs; return false; }

		if (size > rsize)
		{
			g_sub(_d, _d, size, rhs._d, rsize);
			_norm();
			return true;
		}
		if (size < rsize)
		{
			_reallocate(rsize);
			g_sub(_d, rhs._d, rsize, _d, size);
			_norm();
			return false;
		}

		uint64_t * const d = _d;
		const int cmp = g_cmp(d, rhs._d, size);
		if (cmp > 0)
		{
			g_sub(d, d, size, rhs._d, size);
			_norm();
			return true;
		}
		if (cmp < 0)
		{
			g_sub(d, rhs._d, size, d, size);
			_norm();
			return false;
		}
		_size = 0;
		return true;
	}

	// *this >= rhs is required
	guint & operator-=(const guint & rhs) { sub(rhs); return *this; }

	guint & mul(const guint & x, const guint & y)	// *this != x, *this != y
	{
		const size_t xsize = x._size, ysize = y._size;
		if ((xsize == 0) || (ysize == 0)) { *this = 0u; return *this; }
		const size_t new_size = xsize + ysize;
		_reallocate(new_size);
		uint64_t * const d = _d;
		g_mul(d, x._d, xsize, y._d, ysize);
		if (d[new_size - 1] == 0) --_size;
		return *this;
	}

	guint & sqr(const guint & x)	// *this != x
	{
		const size_t size = x._size;
		if (size == 0) { *this = 0u; return *this; }
		const size_t new_size = size + size;
		_reallocate(new_size);
		uint64_t * const d = _d;
		g_sqr(d, x._d, size);
		if (d[new_size - 1] == 0) --_size;
		return *this;
	}

	guint & operator*=(const guint & rhs) { guint t; t.mul(*this, rhs); swap(t); return *this; }

	guint & lshift(const size_t n)
	{
		const size_t size = _size;
		if (size != 0)
		{
			_reallocate(size + n);
			uint64_t * const d = _d;
			g_copy_rev(&d[n], d, size);
			g_zero(d, n);
		}
		return *this;
	}

	guint & rshift(const size_t n)
	{
		const size_t size = _size;
		if (size <= n) *this = 0u;
		else
		{
			g_copy(_d, &_d[n], size - n);
			_reallocate(size - n);
		}
		return *this;
	}

	guint & operator<<=(const size_t s)
	{
		const size_t size = _size;
		if (size == 0) return *this;
		const size_t s_size = s / 64, s_64 = s % 64;
		if (s_64 == 0) return lshift(s_size);

		const uint64_t h = _d[size - 1] >> (64 - s_64);
		_reallocate(size + s_size + ((h != 0) ? 1 : 0));
		uint64_t * const d = _d;
		if (h != 0) d[size + s_size] = h;
		for (size_t j = size - 1; j > 0; --j) d[j + s_size] = (d[j] << s_64) | (d[j - 1] >> (64 - s_64));
		d[s_size] = d[0] << s_64;
		g_zero(d, s_size);
		return *this;
	}

	guint & operator>>=(const size_t s)
	{
		const size_t s_size = s / 64, s_64 = s % 64, size = _size;
		if (size <= s_size) { *this = 0u; return *this; }
		if (s_64 == 0) return rshift(s_size);

		uint64_t * const d = _d;
		for (size_t i = 0; i < size - s_size - 1; ++i) d[i] = (d[i + s_size] >> s_64) | (d[i + s_size + 1] << (64 - s_64));
		const uint64_t l = d[size - 1] >> s_64;
		d[size - s_size - 1] = l;
		_reallocate(size - s_size - ((l == 0) ? 1 : 0));
		return *this;
	}

	guint & div_norm(int & right_shift)
	{
		if (_size == 0) throw std::runtime_error("divide by zero");

		int rshift = 0;
		const uint64_t * const d = _d;
		for (size_t i = 0; d[i] == 0; ++i) rshift += 8 * sizeof(uint64_t);
		uint64_t hi = d[_size - 1]; while (hi >> 63 != 1) { hi *= 2; --rshift; }
		right_shift = rshift;

		if (rshift > 0) *this >>= rshift; else if (rshift < 0) *this <<= -rshift;
		return *this;
	}

	// 2^{2*d_bit_size} / d: d_inv > d
	guint & div_invert(const guint & d)	// *this != d
	{
		const size_t dsize = d._size;
		if (dsize == 0) throw std::runtime_error("divide by zero");

// std::cout << d.to_string() << " (" << dsize << ")" << std::endl;

		guint x;
		x._reallocate(dsize + 1);
		g_zero(x._d, dsize - 1);
		x._d[dsize] = 1; x._d[dsize - 1] = (~size_t(0) / (d._d[dsize - 1] >> 32)) << 32;

// std::cout << x.to_string() << " (" << x._size << ")" << std::endl;

		*this = 0u;
		size_t h = 0;
		while (cmp(x) != 0)
		{
			*this = x;
			guint t; t.sqr(x); t.rshift(dsize - 1); t *= d; t.rshift(dsize + 1);
			x += x; x -= t;
			h++;
			// if (h == 100) break;
		}
		*this -= 1u;

// std::cout << to_string() << " (" << _size << "): " << h << std::endl;
// exit(0);

		return *this;
	}
	// *this: dividend, d: divisor, q: quotient, r: remainder
	// Sources and destinations must be different
	void div_rem(guint & q, guint & r, const guint & d, const guint & d_inv) const
	{
		const size_t size = _size, dsize = d._size;

		if (size < dsize) { q = 0u; r = *this; return; }

		const size_t n = size / dsize + 1;	// >= 2
		r = 0u, q = 0u;
		guint q_j, t;
		for (size_t i = 0, j = n - 1; i < n; ++i, --j)
		{
			if (r._size != 0)
			{
				q_j = r; q_j.lshift(dsize);
				g_copy(q_j._d, &_d[j * dsize], dsize);
			}
			else
			{
				// TODO: first loop, external.
				q_j._reallocate(dsize);
				for (size_t k = 0; k < dsize; ++k) q_j._d[k] = (j * dsize + k < size) ? _d[j * dsize + k] : 0;
				q_j._norm();
			}

			// TODO: single variable x in, r out
			t._div(r, q_j, d, d_inv);

			// TODO: set block
			q.lshift(dsize); q += t;
		}
	}

	// *this /= d, remainder must be zero
	guint & div_exact(const guint & d, const guint & d_inv, const int right_shift)
	{
		if (right_shift > 0) *this >>= right_shift; else if (right_shift < 0) *this <<= -right_shift;

		guint q, r; div_rem(q, r, d, d_inv);
		if (r._size != 0) throw std::runtime_error("div_exact failed");
 		*this = q;
		return *this;
	}

	// *this = x / y, fast if x / y is small
	guint & quotient(const guint & x, const guint & y)	// *this != x, *this != y
	{
		const size_t xsize = x._size, ysize = y._size;

		if (ysize == 0) throw std::runtime_error("divide by zero");
		if ((xsize == 0) || (ysize > xsize)) { *this = 0u; return *this; }

		if (xsize == 1) { *this = x._d[0] / y._d[0]; return *this; }

		if (xsize == ysize)	// Lehmer's gcd: 97.3%
		{
			__uint128_t n = (__uint128_t(x._d[xsize - 1]) << 64) | x._d[xsize - 2];
			__uint128_t d = (__uint128_t(y._d[ysize - 1]) << 64) | y._d[ysize - 2];
			if ((n == ~__uint128_t(0)) || (d == ~__uint128_t(0))) { n >>= 1; d >>= 1; }
			const __uint128_t q_1 = (n + 1) / d, q_2 = n / (d + 1);
			if (q_1 == q_2)
			{
				const uint64_t q_h = uint64_t(q_1 >> 64), q_l = uint64_t(q_1);
				if (q_h != 0) { _reallocate(2); _d[0] = q_l; _d[1] = q_h; }
				else *this = q_l;
				return *this;
			}
		}

		if (xsize == ysize + 1)	// Lehmer's gcd: 2.7%
		{
			__uint128_t n = (__uint128_t(x._d[xsize - 1]) << 64) | x._d[xsize - 2];
			uint64_t d = y._d[ysize - 1];
			if (n == ~__uint128_t(0)) { n >>= 1; d >>= 1; }
			const __uint128_t q_1 = (n + 1) / d, q_2 = n / (__uint128_t(d) + 1);
			if (q_1 == q_2)
			{
				const uint64_t q_h = uint64_t(q_1 >> 64), q_l = uint64_t(q_1);
				if (q_h != 0) { _reallocate(2); _d[0] = q_l; _d[1] = q_h; }
				else *this = q_l;
				return *this;
			}
		}

		// Lehmer's gcd: 0%
		int right_shift; guint d = y, d_inv; d.div_norm(right_shift); d_inv.div_invert(d);
		*this = x; if (right_shift > 0) *this >>= right_shift; else if (right_shift < 0) *this <<= -right_shift;
		guint q, r; div_rem(q, r, d, d_inv);
		*this = q;
		return *this;
	}

	void split(guint & lo, const size_t n)
	{
		const size_t size = _size;
		if (n >= size) throw std::runtime_error("split failed");

		lo._reallocate(n);
		g_copy(lo._d, _d, n);
		lo._norm();

		g_copy(_d, &_d[n], size - n);
		_reallocate(size - n);
	}

	gfloat to_float() const
	{
		const size_t size = _size;
		if (size == 0) return gfloat(0, 0);

		// base 2
		long double mantissa; size_t exponent;
		if (size == 1) { mantissa = _d[0]; exponent = 0; }
		else { mantissa = std::ldexpl(_d[size - 1], 64) + _d[size - 2]; exponent = (size - 2) * sizeof(uint64_t) * 8; };
		while (mantissa >= 1) { mantissa *= 0.5; ++exponent; }
		return gfloat(mantissa, exponent);
	}

	std::string to_string() const
	{
		char * const cstr = new char[_size * 20 + 16];
		g_get_str(cstr, _d, _size);
		const std::string str = std::string(cstr);
		delete[] cstr;
		return str;
	}
};
