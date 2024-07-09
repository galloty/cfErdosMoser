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

// giant integer
class gint
{
private:
	size_t _alloc;
	size_t _size;
	uint64_t * _d;
	bool _is_positive;

private:
	void _allocate(const size_t size)
	{
		_alloc = (size / 1024 + 1) * 1024;
		_size = size;
		_d = (uint64_t *)Heap::allocate_function(_alloc * sizeof(uint64_t));
	}

	void _reallocate(const size_t size, const bool forced = false)
	{
		if ((size > _alloc) || (_alloc - size > 1024 * 1024) || forced)
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

		const int sgn = is_positive ? 1 : -1;
		if (is_positive != ris_positive) return sgn;
		if (size != rsize) return (size > rsize) ? sgn : -sgn;
		const int ucmp = g_cmp(_d, rhs._d, size);
		return is_positive ? ucmp : -ucmp;
	}

	void _uadd(const uint64_t n)
	{
		const uint64_t carry = g_add_1(_d, _size, n);
		if (carry != 0)
		{
			_reallocate(_size + 1);
			_d[_size - 1] = carry;
		}
	}

	void _usub(const uint64_t n)
	{
		const size_t size = _size;
		uint64_t * const d = _d;

		if ((size == 1) && (d[0] < n))	// borrow
		{
			d[0] = n - d[0];
			_is_positive = false;
		}
		else
		{
			g_sub_1(d, size, n);
			if (d[size - 1] == 0) --_size;
		}
	}

	void _umul(const uint64_t n)
	{
		const uint64_t carry = g_mul_1(_d, _size, n);
		if (carry != 0)
		{
			_reallocate(_size + 1);
			_d[_size - 1] = carry;
		}
	}

	void _uaddu(const gint & rhs)	// this and rhs are positive
	{
		const size_t size = _size, rsize = rhs._size;

		_reallocate(std::max(size, rsize) + 1);
		uint64_t * const d = _d;
		const uint64_t carry = (size >= rsize) ? g_add(d, d, size, rhs._d, rsize) : g_add(d, rhs._d, rsize, d, size);
		if (carry != 0) d[_size - 1] = carry; else --_size;
	}

	void _usubu(const gint & rhs)	// this and rhs are positive
	{
		const size_t size = _size, rsize = rhs._size;

		if (size > rsize)
		{
			g_sub(_d, _d, size, rhs._d, rsize);
			_norm();
		}
		else if (size < rsize)
		{
			_reallocate(rsize);
			g_sub(_d, rhs._d, rsize, _d, size);
			_norm();
			_is_positive = false;
		}
		else
		{
			uint64_t * const d = _d;
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

	void _umulu(const gint & x, const gint & y)	// *this != x, *this != y
	{
		const size_t xsize = x._size, ysize = y._size;

		_reallocate(xsize + ysize);
		const uint64_t carry = (xsize >= ysize) ? g_mul(_d, x._d, xsize, y._d, ysize) : g_mul(_d, y._d, ysize, x._d, xsize);
		if (carry == 0) --_size;
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
		_allocate(0);
		_is_positive = true;
	}

	gint(const gint & rhs)
	{
		_allocate(rhs._size);
		_is_positive = rhs._is_positive;
		g_copy(_d, rhs._d, _size);
	}

	gint(const uint64_t n)
	{
		_allocate((n == 0) ? 0 : 1);
		_is_positive = true;
		_d[0] = n;
	}

	virtual ~gint() { _free(); }

	void reset()
	{
		_reallocate(0, true);
		_is_positive = true;
	}

	size_t get_word_count() const { return _size; }
	size_t get_byte_count() const { return get_word_count() * sizeof(uint64_t); }

	gint & operator = (const uint64_t n)
	{
		_reallocate((n == 0) ? 0 : 1);
		_is_positive = true;
		_d[0] = n;
		return *this;
	}

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;

		_reallocate(rhs._size);
		_is_positive = rhs._is_positive;
		g_copy(_d, rhs._d, rhs._size);
		return *this;
	}

	gint & swap(gint & rhs)
	{
		std::swap(_alloc, rhs._alloc);
		std::swap(_size, rhs._size);
		std::swap(_d, rhs._d);
		std::swap(_is_positive, rhs._is_positive);
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

	static void mod_init(const uint64_t n, uint64_t & n_inv, int & n_e) { g_mod_init(n, n_inv, n_e); }

	uint64_t mod(const uint64_t n, const uint64_t n_inv, const int n_e, const uint64_t n_f) const
	{
		if (n == 0) throw std::runtime_error("divide by zero");
		if (_size == 0) return 0;
		const uint64_t remainder = g_mod_1(_d, _size, n, n_inv, n_e, n_f);
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
		_umulu(x, y);
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
		const size_t xsize = x._size, ysize = y._size;

		if (!x._is_positive || !y._is_positive) throw std::runtime_error("divide: negative input");
		if (ysize == 0) throw std::runtime_error("divide by zero");
		if (xsize == 0) { r = 0u; *this = 0u; return *this; }

		if (ysize > xsize) { r = x; *this = 0u; return *this; }

		_reallocate(xsize - ysize + 1); r._reallocate(ysize);
		g_tdiv_qr(_d, r._d, x._d, xsize, y._d, ysize);
		_norm(); r._norm();
		return *this;
	}

	gint & lshift(const size_t n)
	{
		if (_size == 0) return *this;
		_reallocate(_size + n);
		g_copy_rev(&_d[n], _d, _size - n);
		g_zero(_d, n);
		return *this;
	}

	void split(gint & lo, const size_t n)
	{
		if (!_is_positive || (n >= _size)) throw std::runtime_error("split failed");

		lo._reallocate(n);
		g_copy(lo._d, _d, n);
		lo._norm();
		lo._is_positive = true;

		g_copy(_d, &_d[n], _size - n);
		_reallocate(_size - n);
	}

	gfloat to_float() const
	{
		if (_size == 0) return gfloat(0, 0);

		// base 2
		long double mantissa; size_t exponent;
		if (_size == 1) { mantissa = _d[0]; exponent = 0; }
		else { mantissa = std::ldexpl(_d[_size - 1], 64) + _d[_size - 2]; exponent = (_size - 2) * sizeof(uint64_t) * 8; };
		while (mantissa >= 1) { mantissa *= 0.5; ++exponent; }

		return gfloat(_is_positive ? mantissa : -mantissa, exponent);
	}

	std::string to_string() const
	{
		char * const cstr = new char[_size * 20 + 16];
		g_get_str(cstr, _d, _size);
		const std::string str_pos = std::string(cstr);
		delete[] cstr;
		return _is_positive ? str_pos : (std::string("-") + str_pos);
	}
};
