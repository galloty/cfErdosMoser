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

public:
	gint()
	{
		// mpz_init(_z);

		_allocate(1);
		_z->_mp_size = 0;
	}

	virtual ~gint()
	{
		// mpz_clear(_z);

		_free();
	}

	gint(const gint & rhs)
	{
		// mpz_init_set(_z, rhs._z);

		_allocate(std::abs(rhs._z->_mp_size));
		_z->_mp_size = rhs._z->_mp_size;

		for (size_t i = 0, s = std::abs(_z->_mp_size); i < s; ++i) _z->_mp_d[i] = rhs._z->_mp_d[i];
	}

	void reset()
	{
		// mpz_clear(_z); mpz_init(_z);

		_reallocate(1);
		_z->_mp_size = 0;
	}

	size_t get_word_count() const { return std::abs(_z->_mp_size); }
	size_t get_byte_count() const { return get_word_count() * sizeof(mp_limb_t); }

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;
		// mpz_set(_z, rhs._z);

		_reallocate(std::abs(rhs._z->_mp_size));
		_z->_mp_size = rhs._z->_mp_size;

		for (size_t i = 0, s = std::abs(_z->_mp_size); i < s; ++i) _z->_mp_d[i] = rhs._z->_mp_d[i];

		return *this;
	}

	gint & operator = (const uint64_t n)
	{
		// mpz_set_ui_64(_z, n);

		_reallocate(1);
		_z->_mp_size = (n == 0) ? 0 : 1;
		_z->_mp_d[0] = n;
		return *this;
	}

	gint & swap(gint & rhs)
	{
		// mpz_swap(_z, rhs._z);

		std::swap(_z->_mp_alloc, rhs._z->_mp_alloc);
		std::swap(_z->_mp_size, rhs._z->_mp_size);
		std::swap(_z->_mp_d, rhs._z->_mp_d);
		return *this;
	}

	int _cmp(const gint & rhs) const
	{
		if (_z->_mp_size != rhs._z->_mp_size) return (_z->_mp_size > rhs._z->_mp_size) ? 1 : -1;
		const int ucmp = mpn_cmp(_z->_mp_d, rhs._z->_mp_d, std::abs(_z->_mp_size));
		return (std::abs(_z->_mp_size) >= 0 ? ucmp : -ucmp);
	}

	bool operator==(const gint & rhs) const { return (_cmp(rhs) == 0); }
	bool operator!=(const gint & rhs) const { return (_cmp(rhs) != 0); }
	bool operator>=(const gint & rhs) const { return (_cmp(rhs) >= 0); }

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

	gint & operator+=(const uint64_t n)
	{
		// mpz_add_ui(_z, _z, n);

		if (_z->_mp_size == 0) { *this = n; }
		else if (_z->_mp_size > 0) _add(n);
		else { _z->_mp_size = -_z->_mp_size; _sub(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & operator-=(const uint64_t n)
	{
		// mpz_sub_ui(_z, _z, n);

		if (_z->_mp_size == 0) { *this = n; _z->_mp_size = -_z->_mp_size; }
		else if (_z->_mp_size > 0) _sub(n);
		else { _z->_mp_size = -_z->_mp_size; _add(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	gint & operator*=(const uint64_t n)
	{
		// mpz_mul_ui(_z, _z, n);

		if ((_z->_mp_size == 0) || (n == 0)) { *this = 0u; }

		if (_z->_mp_size > 0) _mul(n);
		else { _z->_mp_size = -_z->_mp_size; _mul(n); _z->_mp_size = -_z->_mp_size; }
		return *this;
	}

	uint64_t operator%(const uint64_t n) const
	{
		// return mpz_fdiv_ui(_z, n);

		if (n == 0) throw std::runtime_error("divide by zero");
		if (_z->_mp_size == 0) return 0;
		const mp_limb_t remainder = mpn_mod_1(_z->_mp_d, std::abs(_z->_mp_size), n);
		if (remainder == 0) return 0;
		return (_z->_mp_size < 0) ? n - remainder : remainder;
	}

	gint & operator+=(const gint & rhs) { mpz_add(_z, _z, rhs._z); return *this; }
	gint & operator-=(const gint & rhs) { mpz_sub(_z, _z, rhs._z); return *this; }
	gint & operator*=(const gint & rhs) { mpz_mul(_z, _z, rhs._z); return *this; }
	gint & operator%=(const gint & rhs) { mpz_mod(_z, _z, rhs._z); return *this; }

	// gint & mul(const gint & x, const gint & y) { mpz_mul(_z, x._z, y._z); return *this; }
	gint & addmul(const gint & x, const gint & y) { mpz_addmul(_z, x._z, y._z); return *this; }
	gint & submul(const gint & x, const gint & y) { mpz_submul(_z, x._z, y._z); return *this; }
	gint & div(const gint & x, const gint & y) { mpz_fdiv_q(_z, x._z, y._z); return *this; }

	gint & lshift(const size_t n)
	{
		const size_t size = mpz_size(_z);
		const int sgn = mpz_sgn(_z);
		mp_ptr limbs = mpz_limbs_modify(_z, size + n);
		for (size_t i = 0, j = size - 1; i < size; ++i, --j) limbs[j + n] = limbs[j];
		for (size_t i = 0; i < n; ++i) limbs[i] = 0;
		mpz_limbs_finish(_z, sgn * mp_size_t(size + n));
		return *this;
	}

	void split(gint & lo, const size_t n, const bool fix_roundoff)
	{
		const size_t size = mpz_size(_z);
		if ((mpz_sgn(_z) <= 0) || (n >= size)) throw std::runtime_error("split failed");

		const size_t size_hi = size - n;
		mp_ptr limbs = mpz_limbs_modify(_z, size_hi);

		size_t size_lo = n;
		mp_ptr limbs_lo = mpz_limbs_write(lo._z, size_lo);
		for (size_t i = 0; i < n; ++i) limbs_lo[i] = limbs[i];
		while ((size_lo != 0) && (limbs_lo[size_lo - 1] == 0)) --size_lo;
		mpz_limbs_finish(lo._z, mp_size_t(size_lo));

		for (size_t i = n; i < size; ++i) limbs[i - n] = limbs[i];
		mpz_limbs_finish(_z, mp_size_t(size_hi));

		if (fix_roundoff) mpz_add_ui(_z, _z, 1);
	}

	gint & divexact(const gint & rhs)
	{
		if (mpz_divisible_p(_z, rhs._z)) mpz_divexact(_z, _z, rhs._z); else throw std::runtime_error("divexact failed");
		return *this;
	}

	// uint32_t nu(const uint32_t p) const
	// {
	// 	gint t = *this;
	// 	uint32_t a = 0; while (mpz_divisible_ui_p(t._z, p)) { mpz_divexact_ui(t._z, t._z, p); ++a; }
	// 	return a;
	// }

	gfloat to_float() const
	{
		const size_t size = mpz_size(_z);
		if (size == 0) return gfloat(0, 0);
		const int sgn = mpz_sgn(_z);

		// base 2
		long double mantissa; size_t exponent;
		mp_srcptr limbs = mpz_limbs_read(_z);
		if (size == 1) { mantissa = limbs[0]; exponent = 0; }
		else { mantissa = std::ldexpl(limbs[size - 1], 64) + limbs[size - 2]; exponent = (size - 2) * sizeof(mp_limb_t) * 8; };
		while (mantissa >= 1) { mantissa *= 0.5; ++exponent; }

		return gfloat(sgn * mantissa, exponent);
	}

	std::string to_string() const
	{
		char * const cstr = new char[mpz_sizeinbase(_z, 10) + 16];
		mpz_get_str(cstr, 10, _z);
		const std::string str(cstr);
		delete[] cstr;
		return str;
	}
};
