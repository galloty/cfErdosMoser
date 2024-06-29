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

	static void * allocate_function(size_t size) { _size += size; _max_size = std::max(_max_size, _size); return malloc(size); }
	static void * reallocate_function(void *ptr, size_t old_size, size_t new_size) { _size += new_size - old_size; _max_size = std::max(_max_size, _size); return realloc(ptr, new_size); }
	static void free_function(void *ptr, size_t size) { _size -= size; _max_size = std::max(_max_size, _size); free(ptr); }

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
	mpz_t _z;

public:
	gint() { mpz_init(_z); }
	virtual ~gint() { mpz_clear(_z); }
	gint(const gint & rhs) { mpz_init_set(_z, rhs._z); }

	void reset() { mpz_clear(_z); mpz_init(_z); }	// Free memory

	size_t get_word_count() const { return mpz_size(_z); }
	size_t get_byte_count() const { return get_word_count() * sizeof(mp_limb_t); }

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;
		mpz_set(_z, rhs._z);
		return *this;
	}

	gint & operator = (const uint32_t n) { mpz_set_ui(_z, n); return *this; }
	gint & operator = (const uint64_t n) { mpz_set_ui_64(_z, n); return *this; }

	gint & swap(gint & rhs) { mpz_swap(_z, rhs._z); return *this; }

	int sgn() const { return mpz_sgn(_z); }

	bool operator==(const gint & rhs) const { return (mpz_cmp(_z, rhs._z) == 0); }
	bool operator!=(const gint & rhs) const { return (mpz_cmp(_z, rhs._z) != 0); }
	bool operator>=(const gint & rhs) const { return (mpz_cmp(_z, rhs._z) >= 0); }

	gint & operator+=(const uint32_t n) { mpz_add_ui(_z, _z, n); return *this; }
	gint & operator-=(const uint32_t n) { mpz_sub_ui(_z, _z, n); return *this; }
	gint & operator*=(const uint32_t n) { mpz_mul_ui(_z, _z, n); return *this; }
	uint32_t operator%(const uint32_t n) const { return mpz_fdiv_ui(_z, n); }

	gint & operator+=(const gint & rhs) { mpz_add(_z, _z, rhs._z); return *this; }
	gint & operator*=(const gint & rhs) { mpz_mul(_z, _z, rhs._z); return *this; }

	gint & mul(const gint & x, const gint & y) { mpz_mul(_z, x._z, y._z); return *this; }
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

	// *this must be positive and different from hi & lo
	// if fix_roundoff then hi += 1 and lo -= GMP_LIMB_BITS^n
	void split(gint & hi, gint & lo, const size_t n, const bool fix_roundoff) const
	{
		const size_t size = mpz_size(_z);
		if ((mpz_sgn(_z) <= 0) || (n >= size)) throw std::runtime_error("split failed");

		mp_srcptr limbs = mpz_limbs_read(_z);

		size_t size_lo = n;
		mpz_set_ui(lo._z, 1); mp_ptr limbs_lo = mpz_limbs_write(lo._z, size_lo);
		for (size_t i = 0; i < n; ++i) limbs_lo[i] = limbs[i];
		while ((size_lo != 0) && (limbs_lo[size_lo - 1] == 0)) --size_lo;
		mpz_limbs_finish(lo._z, mp_size_t(size_lo));

		if (fix_roundoff)
		{
			mpz_set_ui(hi._z, 1); mp_ptr limbs_hi = mpz_limbs_write(hi._z, n + 1);
			for (size_t i = 0; i < n; ++i) limbs_hi[i] = 0;
			limbs_hi[n] = 1;
			mpz_limbs_finish(hi._z, mp_size_t(n + 1));
			mpz_sub(lo._z, lo._z, hi._z);
		}

		const size_t size_hi = size - n;
		mpz_set_ui(hi._z, 1); mp_ptr limbs_hi = mpz_limbs_write(hi._z, size_hi);
		for (size_t i = n; i < size; ++i) limbs_hi[i - n] = limbs[i];
		mpz_limbs_finish(hi._z, mp_size_t(size_hi));

		if (fix_roundoff) mpz_add_ui(hi._z, hi._z, 1);
	}

	gint & divexact(const gint & rhs)
	{
		if (mpz_divisible_p(_z, rhs._z)) mpz_divexact(_z, _z, rhs._z); else throw std::runtime_error("divexact failed");
		return *this;
	}

	uint32_t nu(const uint32_t p) const
	{
		gint t = *this;
		uint32_t a = 0; while (mpz_divisible_ui_p(t._z, p)) { mpz_divexact_ui(t._z, t._z, p); ++a; }
		return a;
	}

	void get_d_2exp(double & mantissa, long & exponent) const { mantissa = mpz_get_d_2exp(&exponent, _z); }

	std::string to_string() const
	{
		char * const cstr = new char[mpz_sizeinbase(_z, 10) + 16];
		mpz_get_str(cstr, 10, _z);
		const std::string str(cstr);
		delete[] cstr;
		return str;
	}
};
