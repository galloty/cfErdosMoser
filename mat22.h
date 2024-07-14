/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>

#include "gint.h"

class Mat22
{
private:
	gint _a11, _a12;
	gint _a21, _a22;

public:
	Mat22() {}
	virtual ~Mat22() {}
	Mat22(const Mat22 & rhs) : _a11(rhs._a11), _a12(rhs._a12), _a21(rhs._a21), _a22(rhs._a22) {}

	Mat22 & operator=(const Mat22 & rhs)
	{
		if (&rhs == this) return *this;
		_a11 = rhs._a11; _a12 = rhs._a12; _a21 = rhs._a21; _a22 = rhs._a22;
		return *this;
	}

	const gint & get11() const { return _a11; }
	const gint & get12() const { return _a12; }
	const gint & get21() const { return _a21; }
	const gint & get22() const { return _a22; }

	void set_zero() { _a11 = 0u; _a12 = 0u; _a21 = 0u; _a22 = 0u; }
	void set_identity() { _a11 = 1u; _a12 = 0u; _a21 = 0u; _a22 = 1u; }

	void set_gcf(const uint64_t n)
	{
		_a11 = n; _a12 = 2 * n + 1; _a12 *= _a11;
		_a21 = 2u; _a22 = 5 * n + 2;
	}

	size_t get_min_word_count() const { return std::min(std::min(_a11.get_word_count(), _a12.get_word_count()), std::min(_a21.get_word_count(), _a22.get_word_count())); }
	size_t get_byte_count() const { return _a11.get_byte_count() + _a12.get_byte_count() + _a21.get_byte_count() + _a22.get_byte_count(); }

	Mat22 & swap(Mat22 & rhs)
	{
		_a11.swap(rhs._a11); _a12.swap(rhs._a12);
		_a21.swap(rhs._a21); _a22.swap(rhs._a22);
		return *this;
	}

	Mat22 & mul_right(const Mat22 & rhs)
	{
		gint t1, t2;

		t1.mul(_a11, rhs._a12); t2.mul(_a12, rhs._a22); t1 += t2; t2 = 0u;
		_a11 *= rhs._a11; _a12 *= rhs._a21; _a11 += _a12;
		_a12.swap(t1);

		t1.mul(_a21, rhs._a12); t2.mul(_a22, rhs._a22); t1 += t2; t2 = 0u;
		_a21 *= rhs._a11; _a22 *= rhs._a21; _a21 += _a22;
		_a22.swap(t1);

		return *this;
	}

	Mat22 & div(gint & d)
	{
		int right_shift; d.div_norm(right_shift);
		gint d_inv; d_inv.div_invert(d);
		_a11.div_exact(d, d_inv, right_shift); _a12.div_exact(d, d_inv, right_shift);
		_a21.div_exact(d, d_inv, right_shift); _a22.div_exact(d, d_inv, right_shift);
		return *this;
	}

	Mat22 & mul_left(const Mat22 & rhs)
	{
		gint t1, t2;

		t1.mul(_a11, rhs._a21); t2.mul(_a21, rhs._a22); t1 += t2;
		_a11 *= rhs._a11; _a21 *= rhs._a12; _a11 += _a21;
		_a21.swap(t1);

		t1.mul(_a12, rhs._a21); t2.mul(_a22, rhs._a22); t1 += t2;
		_a12 *= rhs._a11; _a22 *= rhs._a12; _a12 += _a22;
		_a22.swap(t1);

		return *this;
	}

	void split(Mat22 & lo, const size_t n)
	{
		// We may have hi.a_11/hi.a_21 > a_11/a_21 and hi.a_12/hi.a_22 < a_12/a_22 if hi.a_i = a_i >> (n * GMP_LIMB_BITS)
		// If hi.a_21 = (a_21 >> (n * GMP_LIMB_BITS)) + 1 then hi.a_11/hi.a_21 < a_11/a_21.
		// If hi.a_12 = (a_12 >> (n * GMP_LIMB_BITS)) + 1 then hi.a_12/hi.a_22 > a_12/a_22.
		// Since a_11/a_21 < alpha < a_12/a_22, we still have hi.a_11/hi.a_21 < alpha < hi.a_12/hi.a_22.
		_a11.split(lo._a11, n); _a12.split(lo._a12, n); _a12 += 1u;
		_a21.split(lo._a21, n); _a21 += 1u; _a22.split(lo._a22, n);
	}

	void combine(const Mat22 & lo, const Mat22 & Ml, const size_t n)
	{
		// M -= Ml * [0 1 1 0]
		_a11 -= Ml._a12; _a11.lshift(n); _a11 += lo._a11;
		_a12 -= Ml._a11; _a12.lshift(n); _a12 += lo._a12;
		_a21 -= Ml._a22; _a21.lshift(n); _a21 += lo._a21;
		_a22 -= Ml._a21; _a22.lshift(n); _a22 += lo._a22;
	}

	void init_gcf(const gint & N)
	{
		_a11 = 0u; _a12 = 1u;
		_a21 = N; _a21 *= 2u; _a22 = _a21;
	}

	void cf_mul(const gint & coefficient)
	{
		gint t;
		t.mul(coefficient, _a21); _a11 -= t; _a11.swap(_a21);
		t.mul(coefficient, _a22); _a12 -= t; _a12.swap(_a22);
	}

	void cf_mul_revert(const gint & coefficient)
	{
		gint t;
		_a11.swap(_a21); t.mul(coefficient, _a21); _a11 += t;
		_a12.swap(_a22); t.mul(coefficient, _a22); _a12 += t;
	}

	bool get_cf_coefficient(gint & coefficient1, gint & coefficient2)
	{
		gint t;

		coefficient1.quotient(_a11, _a21); t.quotient(_a12, _a22);
		if (coefficient1 != t) return false;
		cf_mul(coefficient1);

		coefficient2.quotient(_a11, _a21); t.quotient(_a12, _a22);
		const bool success = (coefficient2 == t);
		if (success) cf_mul(coefficient2); else cf_mul_revert(coefficient1);
		return success;
	}
};
