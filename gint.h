/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

#include <gmp.h>

#include "gfloat.h"

// Low-level implemented using GMP functions. But on Windows, mpn has limit of 2^(31 + 6) bits (41 billion digits).

inline int g_cmp(const mp_limb_t * const s1p, const mp_limb_t * const s2p, const size_t n)
{
	return mpn_cmp(s1p, s2p, mp_size_t(n));
}

inline mp_limb_t g_add_1(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t n, const mp_limb_t s2limb)
{
	return mpn_add_1(rp, s1p, mp_size_t(n), s2limb);
}

inline mp_limb_t g_sub_1(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t n, const mp_limb_t s2limb)
{
	return mpn_sub_1(rp, s1p, mp_size_t(n), s2limb);
}

inline mp_limb_t g_mul_1(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t n, const mp_limb_t s2limb)
{
	return mpn_mul_1(rp, s1p, mp_size_t(n), s2limb);
}

inline mp_limb_t g_mod_1(const mp_limb_t * const s1p, const mp_size_t s1n, const mp_limb_t s2limb)
{
	return mpn_mod_1(s1p, mp_size_t(s1n), s2limb);
}

inline mp_limb_t g_add(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t s1n, const mp_limb_t * const s2p, const size_t s2n)
{
	return mpn_add(rp, s1p, mp_size_t(s1n), s2p, mp_size_t(s2n));
}

inline mp_limb_t g_sub(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t s1n, const mp_limb_t * const s2p, const size_t s2n)
{
	return mpn_sub(rp, s1p, mp_size_t(s1n), s2p, mp_size_t(s2n));
}

inline mp_limb_t g_mul(mp_limb_t * const rp, const mp_limb_t * const s1p, const size_t s1n, const mp_limb_t * const s2p, const size_t s2n)
{
	return mpn_mul(rp, s1p, mp_size_t(s1n), s2p, mp_size_t(s2n));
}

inline void g_tdiv_qr(mp_limb_t * const qp, mp_limb_t * const rp, const size_t qxn, const mp_limb_t * const np, const size_t nn, const mp_limb_t * const dp, const size_t dn)
{
	mpn_tdiv_qr(qp, rp, mp_size_t(qxn), np, mp_size_t(nn), dp, mp_size_t(dn));
}

inline size_t g_get_str(unsigned char * const str, mp_limb_t * const s1p, const size_t s1n)
{
	return (size_t)mpn_get_str(str, 10, s1p, mp_size_t(s1n));
}

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

// giant integer
class gint
{
private:
	size_t _alloc;
	size_t _size;
	bool _is_positive;
	mp_limb_t * _d;

private:
	void _allocate(const size_t size)
	{
		_alloc = (size / 1024 + 1) * 1024;
		_d = (mp_limb_t *)Heap::allocate_function(_alloc * sizeof(mp_limb_t));
	}

	void _reallocate(const size_t size, const bool forced = false)
	{
		if ((size > _alloc) || (_alloc - size > 1024 * 1024) || forced)
		{
			const size_t alloc = (size / 1024 + 1) * 1024;
			_d = (mp_limb_t *)Heap::reallocate_function(_d, _alloc * sizeof(mp_limb_t), alloc * sizeof(mp_limb_t));
			_alloc = alloc;
		}
	}

	void _free()
	{
		Heap::free_function(_d, _alloc * sizeof(mp_limb_t));
	}

	void _norm()
	{
		const mp_limb_t * const d = _d;
		size_t size = _size; while ((size > 0) && (d[size - 1] == 0)) --size;
		_size = size;
	}

	static void _copy(mp_limb_t * const dst, const mp_limb_t * const src, const size_t size)
	{
		for (size_t i = 0; i < size; ++i) dst[i] = src[i];
	}

	void _neg() { _is_positive = !_is_positive; }

	int _cmp(const gint & rhs) const
	{
		const size_t size = _size, rsize = rhs._size;
		const bool is_positive = _is_positive, ris_positive = rhs._is_positive;

		// 0 may be positive or negative
		if (size == 0)
		{
			if (rsize == 0) return 0;
			return ris_positive ? -1 : 1;
		}
		if (rsize == 0) return is_positive ? 1 : -1;

		if (is_positive != ris_positive) return is_positive ? 1 : -1;
		const int sgn = is_positive ? 1 : -1;
		if (size != rsize) return (size > rsize) ? sgn : -sgn;
		const int ucmp = g_cmp(_d, rhs._d, size);
		return is_positive ? ucmp : -ucmp;
	}

	void _uadd(const uint64_t n)
	{
		const size_t size = _size;
		mp_limb_t * const d = _d;

		const mp_limb_t carry = g_add_1(d, d, size, n);
		if (carry != 0)
		{
			_reallocate(size + 1);
			_d[size] = carry;
			++_size;
		}
	}

	void _usub(const uint64_t n)
	{
		const size_t size = _size;
		mp_limb_t * const d = _d;

		if ((size == 1) && (d[0] < n))	// borrow
		{
			d[0] = n - d[0];
			_is_positive = false;
		}
		else
		{
			g_sub_1(d, d, size, n);
			if (d[size - 1] == 0) --_size;
		}
	}

	void _umul(const uint64_t n)
	{
		const size_t size = _size;
		mp_limb_t * const d = _d;

		const mp_limb_t carry = g_mul_1(d, d, size, n);
		if (carry != 0)
		{
			_reallocate(size + 1);
			_d[size] = carry;
			++_size;
		}
	}

	void _uaddu(const gint & rhs)	// this and rhs are positive
	{
		const size_t size = _size, rsize = rhs._size;

		const size_t new_size = std::max(size, rsize);
		_reallocate(new_size + 1);

		mp_limb_t * const d = _d;
		const mp_limb_t carry = (size >= rsize) ? g_add(d, d, size, rhs._d, rsize) : g_add(d, rhs._d, rsize, d, size);
		_size = new_size;
		if (carry != 0)
		{
			d[size] = carry;
			++_size;
		}
	}

	void _usubu(const gint & rhs)	// this and rhs are positive
	{
		const size_t size = _size, rsize = rhs._size;

		if (size > rsize)
		{
			mp_limb_t * const d = _d;
			g_sub(d, d, size, rhs._d, rsize);
			if (d[size - 1] == 0) --_size;
		}
		else if (size < rsize)
		{
			_reallocate(rsize);
			mp_limb_t * const d = _d;
			g_sub(d, rhs._d, rsize, d, size);
			_size = rsize;
			if (d[rsize - 1] == 0) _size -= 1;
			_is_positive = false;
		}
		else
		{
			mp_limb_t * const d = _d;
			const int cmp = g_cmp(d, rhs._d, size);
			if (cmp > 0)
			{
				g_sub(d, d, size, rhs._d, size);
				_norm();
			}
			else if (cmp < 0)
			{
				g_sub(d, rhs._d, size, d, size);
				_norm();
				_is_positive = false;
			}
			else _size = 0;
		}
	}

	void _uadd(const gint & rhs)
	{
		if (rhs._size == 0) return;
		if (rhs._is_positive) _uaddu(rhs); else _usubu(rhs);
	}

	void _usub(const gint & rhs)
	{
		if (rhs._size == 0) return;
		if (rhs._is_positive) _usubu(rhs); else _uaddu(rhs);
	}

public:
	gint()
	{
		_allocate(1);
		_size = 0;
		_is_positive = true;
	}

	virtual ~gint()
	{
		_free();
	}

	gint(const gint & rhs)
	{
		const size_t size = rhs._size;
		_allocate(std::max(size, size_t(1)));
		_size = size;
		_is_positive = rhs._is_positive;

		_copy(_d, rhs._d, size);
	}

	gint(const uint64_t n)
	{
		_allocate(1);
		_size = (n == 0) ? 0 : 1;
		_is_positive = true;
		_d[0] = n;
	}

	void reset()
	{
		_reallocate(1, true);
		_size = 0;
		_is_positive = true;
	}

	size_t get_word_count() const { return _size; }
	size_t get_byte_count() const { return get_word_count() * sizeof(mp_limb_t); }

	gint & operator = (const uint64_t n)
	{
		_reallocate(1);
		_size = (n == 0) ? 0 : 1;
		_is_positive = true;
		_d[0] = n;
		return *this;
	}

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;

		const size_t size = rhs._size;
		_reallocate(std::max(size, size_t(1)));
		_size = size;
		_is_positive = rhs._is_positive;

		_copy(_d, rhs._d, size);

		return *this;
	}

	gint & swap(gint & rhs)
	{
		std::swap(_alloc, rhs._alloc);
		std::swap(_size, rhs._size);
		std::swap(_is_positive, rhs._is_positive);
		std::swap(_d, rhs._d);
		return *this;
	}

	bool is_zero() const { return (_size == 0); }
	bool operator==(const gint & rhs) const { return (_cmp(rhs) == 0); }
	bool operator!=(const gint & rhs) const { return (_cmp(rhs) != 0); }
	bool operator>=(const gint & rhs) const { return (_cmp(rhs) >= 0); }

	gint & operator+=(const uint64_t n)
	{
		if (_size == 0) { *this = n; }
		else if (_is_positive) _uadd(n);
		else { _neg(); _usub(n); _neg(); }
		return *this;
	}

	gint & operator-=(const uint64_t n)
	{
		if (_size == 0) { *this = n; _is_positive = false; }
		else if (_is_positive) _usub(n); else _uadd(n);
		return *this;
	}

	gint & operator*=(const uint64_t n)
	{
		if ((_size == 0) || (n == 0)) { *this = 0u; }
		_umul(n);
		return *this;
	}

	uint64_t operator%(const uint64_t n) const
	{
		if (n == 0) throw std::runtime_error("divide by zero");
		if (_size == 0) return 0;
		const mp_limb_t remainder = g_mod_1(_d, _size, n);
		if (remainder == 0) return 0;
		return _is_positive ? remainder : n - remainder;
	}

	gint & operator+=(const gint & rhs)
	{
		if (_size == 0) { *this = rhs; }
		else if (_is_positive) _uadd(rhs);
		else { _neg(); _usub(rhs); _neg(); }
		return *this;
	}

	gint & operator-=(const gint & rhs)
	{
		if (_size == 0) { *this = rhs; _neg(); }
		else if (_is_positive) _usub(rhs);
		else { _neg(); _uadd(rhs); _neg(); }
		return *this;
	}

	gint & mul(const gint & x, const gint & y)	// *this != x, *this != y
	{
		if ((x._size == 0) || (y._size == 0)) { *this = 0u; return *this; }

		const size_t new_size = x._size + y._size;
		_reallocate(new_size);
		const mp_limb_t carry = (x._size >= y._size) ? g_mul(_d, x._d, x._size, y._d, y._size) : g_mul(_d, y._d, y._size, x._d, x._size);
		_size = (carry == 0) ? new_size - 1 : new_size;
		_is_positive = (x._is_positive == y._is_positive);
		return *this;
	}

	gint & operator*=(const gint & rhs)
	{
		gint t; t.mul(*this, rhs); swap(t);
		return *this;
	}

	gint & divrem(const gint & x, const gint & y, gint & r)
	{
		if (!x._is_positive || !y._is_positive) throw std::runtime_error("divide: negative input");
		if (y._size == 0) throw std::runtime_error("divide by zero");
		if (x._size == 0) { r = 0u; *this = 0u; return *this; }

		if (y._size > x._size) { r = x; *this = 0u; return *this; }

		const size_t q_size = x._size - y._size + 1, r_size = y._size;

		_reallocate(q_size); r._reallocate(r_size);
		g_tdiv_qr(_d, r._d, 0u, x._d, x._size, y._d, y._size);

		_size = q_size;
		while ((_size > 0) && (_d[_size - 1] == 0)) --_size;

		r._size = r_size;
		while ((r._size > 0) && (r._d[r._size - 1] == 0)) --r._size;

		return *this;
	}

	gint & lshift(const size_t n)
	{
		if (_size == 0) return *this;
		_reallocate(_size + n);
		for (size_t i = 0, j = _size - 1; i < _size; ++i, --j) _d[j + n] = _d[j];
		for (size_t i = 0; i < n; ++i) _d[i] = 0;
		_size += n;
		return *this;
	}

	void split(gint & lo, const size_t n, const bool fix_roundoff)
	{
		if (!_is_positive || (n >= _size)) throw std::runtime_error("split failed");

		size_t size_lo = n;
		lo._reallocate(size_lo);
		for (size_t i = 0; i < n; ++i) lo._d[i] = _d[i];
		while ((size_lo > 0) && (lo._d[size_lo - 1] == 0)) --size_lo;
		lo._size = size_lo;
		lo._is_positive = true;

		const size_t size_hi = _size - n;
		for (size_t i = n; i < _size; ++i) _d[i - n] = _d[i];
		_size = size_hi;
		_reallocate(size_hi);

		if (fix_roundoff) *this += 1u;
	}

	gfloat to_float() const
	{
		if (_size == 0) return gfloat(0, 0);

		// base 2
		long double mantissa; size_t exponent;
		if (_size == 1) { mantissa = _d[0]; exponent = 0; }
		else { mantissa = std::ldexpl(_d[_size - 1], 64) + _d[_size - 2]; exponent = (_size - 2) * sizeof(mp_limb_t) * 8; };
		while (mantissa >= 1) { mantissa *= 0.5; ++exponent; }

		return gfloat(_is_positive ? mantissa : -mantissa, exponent);
	}

	std::string to_string() const
	{
		char * const cstr = new char[_size * 20 + 16];
		const size_t str_size = g_get_str((unsigned char *)cstr, _d, _size);
		for (size_t i = 0; i < str_size; ++i) cstr[i] += '0';
  		cstr[str_size] = '\0';
		const std::string str_pos = std::string(cstr);
		delete[] cstr;
		return _is_positive ? str_pos : (std::string("-") + str_pos);
	}
};
