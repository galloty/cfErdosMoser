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

inline void mpz_set_ui_64(mpz_t rop, const uint64_t n)
{
#ifdef _LONG_LONG_LIMB
	mpz_set_ui(rop, (unsigned long int)n); if (n > 0) rop->_mp_d[0] = n;
#else
	mpz_set_ui(rop, (unsigned long int)(n >> 32));
	mpz_mul_2exp(rop, rop, 32);
	mpz_add_ui(rop, rop, (unsigned long int)n);
#endif
}

class Heap
{
private:
	static size_t _size, _max_size;

public:
	static void * allocate_function(size_t size) { _size += size; _max_size = std::max(_max_size, _size); return malloc(size); }
	static void * reallocate_function(void *ptr, size_t old_size, size_t new_size) { _size += new_size - old_size; _max_size = std::max(_max_size, _size); return realloc(ptr, new_size); }
	static void free_function(void *ptr, size_t size) { _size -= size; _max_size = std::max(_max_size, _size); if (size > 0) free(ptr); }

public:
	Heap() { mp_set_memory_functions(allocate_function, reallocate_function, free_function); }
	virtual ~Heap() { mp_set_memory_functions(nullptr, nullptr, nullptr); }

	size_t get_size() const { return _size; }
	size_t get_max_size() const { return _max_size; }

	void reset_max_size() { _max_size = 0; }
};

size_t Heap::_size = 0, Heap::_max_size = 0;

// giant integer
class gint
{
private:
	// mpz_t has limit of 2^{31 + 6) bits (41 billion digits)
	mpz_t _z;

private:
	void _allocate(const size_t size)
	{
		_z->_mp_alloc = size;
		_z->_mp_d = (mp_limb_t *)Heap::allocate_function(size * sizeof(mp_limb_t));
	}

	void _reallocate(const size_t size)
	{
		// TODO: add block size and logic
		_z->_mp_d = (mp_limb_t *)Heap::reallocate_function(_z->_mp_d, _z->_mp_alloc * sizeof(mp_limb_t), size * sizeof(mp_limb_t));
		_z->_mp_alloc = size;
	}

	void _free()
	{
		Heap::free_function(_z->_mp_d, _z->_mp_alloc * sizeof(mp_limb_t));
	}

	int _cmp(const gint & rhs) const
	{
		if (_z->_mp_size != rhs._z->_mp_size) return (_z->_mp_size > rhs._z->_mp_size) ? 1 : -1;
		const int ucmp = mpn_cmp(_z->_mp_d, rhs._z->_mp_d, std::abs(_z->_mp_size));
		return (std::abs(_z->_mp_size) >= 0 ? ucmp : -ucmp);
	}

	void _add(const uint64_t n)	// _z->_mp_size > 0
	{
		const mp_limb_t carry = mpn_add_1(_z->_mp_d, _z->_mp_d, _z->_mp_size, n);
		if (carry != 0)
		{
			_reallocate(_z->_mp_size + 1);
			_z->_mp_d[_z->_mp_size] = carry;
			_z->_mp_size += 1;
		}
	}

	void _sub(const uint64_t n)	// _z->_mp_size > 0
	{
		if ((_z->_mp_size == 1) && (_z->_mp_d[0] < n))	// borrow
		{
			_z->_mp_d[0] = n - _z->_mp_d[0];
			_z->_mp_size = -1;
		}
		else
		{
			mpn_sub_1(_z->_mp_d, _z->_mp_d, _z->_mp_size, n);
			if (_z->_mp_d[_z->_mp_size - 1] == 0) _z->_mp_size -= 1;
		}
	}

	void _mul(const uint64_t n)	// _z->_mp_size > 0
	{
		const mp_limb_t carry = mpn_mul_1(_z->_mp_d, _z->_mp_d, _z->_mp_size, n);
		if (carry != 0)
		{
			_reallocate(_z->_mp_size + 1);
			_z->_mp_d[_z->_mp_size] = carry;
			_z->_mp_size += 1;
		}
	}

	void _add_pos(const gint & rhs)	// _z->_mp_size > 0 and rhs->_mp_size > 0
	{
		const size_t dst_size = std::max(_z->_mp_size, rhs._z->_mp_size);
		_reallocate(dst_size + 1);
		mp_limb_t carry;
		if (_z->_mp_size >= rhs._z->_mp_size) carry = mpn_add(_z->_mp_d, _z->_mp_d, _z->_mp_size, rhs._z->_mp_d, rhs._z->_mp_size);
		else carry = mpn_add(_z->_mp_d, rhs._z->_mp_d, rhs._z->_mp_size, _z->_mp_d, _z->_mp_size);
		_z->_mp_size = dst_size;
		if (carry != 0)
		{
			_z->_mp_d[_z->_mp_size] = carry;
			_z->_mp_size += 1;
		}
	}

	void _sub_pos(const gint & rhs)	// _z->_mp_size > 0 and rhs->_mp_size > 0
	{
		if (_z->_mp_size > rhs._z->_mp_size)
		{
			mpn_sub(_z->_mp_d, _z->_mp_d, _z->_mp_size, rhs._z->_mp_d, rhs._z->_mp_size);
			if (_z->_mp_d[_z->_mp_size - 1] == 0) _z->_mp_size -= 1;
		}
		else if (_z->_mp_size < rhs._z->_mp_size)
		{
			_reallocate(rhs._z->_mp_size);
			mpn_sub(_z->_mp_d, rhs._z->_mp_d, rhs._z->_mp_size, _z->_mp_d, _z->_mp_size);
			_z->_mp_size = rhs._z->_mp_size;
			if (_z->_mp_d[_z->_mp_size - 1] == 0) _z->_mp_size -= 1;
			_z->_mp_size = -_z->_mp_size;
		}
		else
		{
			const int cmp = mpn_cmp(_z->_mp_d, rhs._z->_mp_d, _z->_mp_size);
			if (cmp > 0)
			{
				mpn_sub_n(_z->_mp_d, _z->_mp_d, rhs._z->_mp_d, _z->_mp_size);
				while ((_z->_mp_size > 0) && (_z->_mp_d[_z->_mp_size - 1] == 0)) _z->_mp_size -= 1;
			}
			else if (cmp < 0)
			{
				mpn_sub_n(_z->_mp_d, rhs._z->_mp_d, _z->_mp_d, _z->_mp_size);
				while ((_z->_mp_size > 0) && (_z->_mp_d[_z->_mp_size - 1] == 0)) _z->_mp_size -= 1;
				_z->_mp_size = -_z->_mp_size;
			}
			else _z->_mp_size = 0;
		}
	}

	void _add(const gint & rhs)	// _z->_mp_size > 0
	{
		if (rhs._z->_mp_size > 0) _add_pos(rhs);
		else if (rhs._z->_mp_size < 0)
		{
			gint t; t._free();
			t._z->_mp_alloc = rhs._z->_mp_alloc; t._z->_mp_d = rhs._z->_mp_d; t._z->_mp_size = -rhs._z->_mp_size;
			_sub_pos(t);
			t._allocate(1);
		}
	}

	void _sub(const gint & rhs)	// _z->_mp_size > 0
	{
		if (rhs._z->_mp_size > 0) _sub_pos(rhs);
		else if (rhs._z->_mp_size < 0)
		{
			gint t; t._free();
			t._z->_mp_alloc = rhs._z->_mp_alloc; t._z->_mp_d = rhs._z->_mp_d; t._z->_mp_size = -rhs._z->_mp_size;
			_add_pos(t);
			t._allocate(1);
		}
	}

public:
	gint()
	{
		_allocate(1);
		_z->_mp_size = 0;
	}

	virtual ~gint()
	{
		_free();
	}

	gint(const gint & rhs)
	{
		_allocate(std::abs(rhs._z->_mp_size));
		_z->_mp_size = rhs._z->_mp_size;

		for (size_t i = 0, s = std::abs(_z->_mp_size); i < s; ++i) _z->_mp_d[i] = rhs._z->_mp_d[i];
	}

	gint(const uint64_t n)
	{
		_allocate(1);
		_z->_mp_size = (n == 0) ? 0 : 1;
		_z->_mp_d[0] = n;
	}

	void reset()
	{
		_reallocate(1);
		_z->_mp_size = 0;
	}

	size_t get_word_count() const { return std::abs(_z->_mp_size); }
	size_t get_byte_count() const { return get_word_count() * sizeof(mp_limb_t); }

	gint & operator = (const uint64_t n)
	{
		_reallocate(1);
		_z->_mp_size = (n == 0) ? 0 : 1;
		_z->_mp_d[0] = n;
		return *this;
	}

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;

		_reallocate(std::abs(rhs._z->_mp_size));
		_z->_mp_size = rhs._z->_mp_size;

		for (size_t i = 0, s = std::abs(_z->_mp_size); i < s; ++i) _z->_mp_d[i] = rhs._z->_mp_d[i];

		return *this;
	}

	gint & swap(gint & rhs)
	{
		std::swap(_z->_mp_alloc, rhs._z->_mp_alloc);
		std::swap(_z->_mp_size, rhs._z->_mp_size);
		std::swap(_z->_mp_d, rhs._z->_mp_d);
		return *this;
	}

	bool is_zero() const { return (_z->_mp_size == 0); }
	bool operator==(const gint & rhs) const { return (_cmp(rhs) == 0); }
	bool operator!=(const gint & rhs) const { return (_cmp(rhs) != 0); }
	bool operator>=(const gint & rhs) const { return (_cmp(rhs) >= 0); }

	gint & operator+=(const uint64_t n)
	{
		if (_z->_mp_size == 0) { *this = n; }
		else if (_z->_mp_size > 0) _add(n);
		else { _z->_mp_size = -_z->_mp_size; _sub(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & operator-=(const uint64_t n)
	{
		if (_z->_mp_size == 0) { *this = n; _z->_mp_size = -_z->_mp_size; }
		else if (_z->_mp_size > 0) _sub(n);
		else { _z->_mp_size = -_z->_mp_size; _add(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & operator*=(const uint64_t n)
	{
		if ((_z->_mp_size == 0) || (n == 0)) { *this = 0u; }

		if (_z->_mp_size > 0) _mul(n);
		else { _z->_mp_size = -_z->_mp_size; _mul(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	uint64_t operator%(const uint64_t n) const
	{
		if (n == 0) throw std::runtime_error("divide by zero");
		if (_z->_mp_size == 0) return 0;
		const mp_limb_t remainder = mpn_mod_1(_z->_mp_d, std::abs(_z->_mp_size), n);
		if (remainder == 0) return 0;
		return (_z->_mp_size < 0) ? n - remainder : remainder;
	}

	gint & operator+=(const gint & rhs)
	{
		if (_z->_mp_size == 0) { *this = rhs; }
		else if (_z->_mp_size > 0) _add(rhs);
		else { _z->_mp_size = -_z->_mp_size; _sub(rhs); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & operator-=(const gint & rhs)
	{
		if (_z->_mp_size == 0) { *this = rhs; _z->_mp_size = -_z->_mp_size; }
		else if (_z->_mp_size > 0) _sub(rhs);
		else { _z->_mp_size = -_z->_mp_size; _add(rhs); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & mul(const gint & x, const gint & y)	// *this != x, *this != y
	{
		if ((x._z->_mp_size == 0) || (y._z->_mp_size == 0)) { *this = 0u; return *this; }

		const size_t x_size = std::abs(x._z->_mp_size), y_size = std::abs(y._z->_mp_size), dst_size = x_size + y_size;
		const int dst_sgn = x._z->_mp_size ^ y._z->_mp_size;
		_reallocate(dst_size);
		mp_limb_t carry;
		if (x_size >= y_size) carry = mpn_mul(_z->_mp_d, x._z->_mp_d, x_size, y._z->_mp_d, y_size);
		else carry = mpn_mul(_z->_mp_d, y._z->_mp_d, y_size, x._z->_mp_d, x_size);

		_z->_mp_size = (carry == 0) ? dst_size - 1 : dst_size;
		if (dst_sgn < 0) _z->_mp_size = -_z->_mp_size;
		return *this;
	}

	gint & operator*=(const gint & rhs)
	{
		gint t; t.mul(*this, rhs); swap(t);
		return *this;
	}

	gint & divrem(const gint & x, const gint & y, gint & r)
	{
		if ((x._z->_mp_size < 0) || (x._z->_mp_size < 0)) throw std::runtime_error("divide: negative input");
		if (y._z->_mp_size == 0) throw std::runtime_error("divide by zero");
		if (x._z->_mp_size == 0) { *this = 0u; return *this; }

		const size_t x_size = std::abs(x._z->_mp_size), y_size = std::abs(y._z->_mp_size);
		if (y_size > x_size) { *this = 0u; return *this; }

		const size_t q_size = x_size - y_size + 1, r_size = y_size;

		_reallocate(q_size); r._reallocate(r_size);
		mpn_tdiv_qr(_z->_mp_d, r._z->_mp_d, 0u, x._z->_mp_d, x_size, y._z->_mp_d, y_size);

		_z->_mp_size = q_size;
		while ((_z->_mp_size > 0) && (_z->_mp_d[_z->_mp_size - 1] == 0)) _z->_mp_size -= 1;

		r._z->_mp_size = r_size;
		while ((r._z->_mp_size > 0) && (r._z->_mp_d[r._z->_mp_size - 1] == 0)) r._z->_mp_size -= 1;

		return *this;
	}

	gint & lshift(const size_t n)
	{
		if (_z->_mp_size == 0) return *this;
		const size_t size = std::abs(_z->_mp_size);
		const int sgn = (_z->_mp_size < 0) ? -1 : 1;
		_reallocate(size + n);
		for (size_t i = 0, j = size - 1; i < size; ++i, --j) _z->_mp_d[j + n] = _z->_mp_d[j];
		for (size_t i = 0; i < n; ++i) _z->_mp_d[i] = 0;
		_z->_mp_size = sgn * mp_size_t(size + n);
		return *this;
	}

	void split(gint & lo, const size_t n, const bool fix_roundoff)
	{
		const size_t size = std::abs(_z->_mp_size);
		if ((_z->_mp_size <= 0) || (n >= size)) throw std::runtime_error("split failed");

		size_t size_lo = n;
		lo._reallocate(size_lo);
		for (size_t i = 0; i < n; ++i) lo._z->_mp_d[i] = _z->_mp_d[i];
		while ((size_lo > 0) && (lo._z->_mp_d[size_lo - 1] == 0)) --size_lo;
		lo._z->_mp_size = size_lo;

		const size_t size_hi = size - n;
		for (size_t i = n; i < size; ++i) _z->_mp_d[i - n] = _z->_mp_d[i];
		_z->_mp_size = size_hi;
		_reallocate(size_hi);

		if (fix_roundoff) *this += 1u;
	}

	gfloat to_float() const
	{
		const size_t size = std::abs(_z->_mp_size);
		if (size == 0) return gfloat(0, 0);
		const int sgn = (_z->_mp_size < 0) ? -1 : 1;

		// base 2
		long double mantissa; size_t exponent;
		if (size == 1) { mantissa = _z->_mp_d[0]; exponent = 0; }
		else { mantissa = std::ldexpl(_z->_mp_d[size - 1], 64) + _z->_mp_d[size - 2]; exponent = (size - 2) * sizeof(mp_limb_t) * 8; };
		while (mantissa >= 1) { mantissa *= 0.5; ++exponent; }

		return gfloat(sgn * mantissa, exponent);
	}

	std::string to_string() const
	{
		char * const cstr = new char[mpn_sizeinbase(_z->_mp_d, _z->_mp_size, 10) + 16];
		const size_t str_size = mpn_get_str((unsigned char *)cstr, 10, _z->_mp_d, _z->_mp_size);
		for (size_t i = 0; i < str_size; ++i) cstr[i] += '0';
  		cstr[str_size] = '\0';
		const std::string str(cstr);
		delete[] cstr;
		return str;
	}
};
