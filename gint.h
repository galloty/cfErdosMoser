/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <iostream>

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

// giant integer
class gint
{
private:
	mpz_t _z;

public:
	gint() { mpz_init(_z); }
	virtual ~gint() { mpz_clear(_z); }
	gint(const gint & rhs) { mpz_init_set(_z, rhs._z); }

	size_t get_word_count() const { return mpz_size(_z) * sizeof(mp_limb_t); }
	size_t get_byte_count() const { return get_word_count() * sizeof(mp_limb_t); }

	gint & operator = (const gint & rhs)
	{
		if (&rhs == this) return *this;
		mpz_set(_z, rhs._z);
		return *this;
	}

	gint & operator = (const uint32_t n) { mpz_set_ui(_z, n); return *this; }
	gint & operator = (const uint64_t n) { mpz_set_ui_64(_z, n); return *this; }

	void swap(gint & rhs) { mpz_swap(_z, rhs._z); }

	bool operator==(const gint & rhs) const { return (mpz_cmp(_z, rhs._z) == 0); }
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

	gint & rshift(const gint & rhs, const size_t s) { mpz_div_2exp(_z, rhs._z, s * 8); return *this; }	// TODO fix: s * 8 may be > 32-bit

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

	void out(FILE * const fp) const { mpz_out_str(fp, 10, _z); }
};
