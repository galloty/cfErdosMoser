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

	Mat22 & operator = (const Mat22 & rhs)
	{
		if (&rhs == this) return *this;
		_a11 = rhs._a11; _a12 = rhs._a12; _a21 = rhs._a21; _a22 = rhs._a22;
		return *this;
	}

	const gint & get11() const { return _a11; }
	const gint & get12() const { return _a12; }
	const gint & get21() const { return _a21; }
	const gint & get22() const { return _a22; }

	void clear() { _a11 = 0u; _a12 = 0u; _a21 = 0u; _a22 = 0u; }
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

	Mat22 & operator+=(const Mat22 & rhs)
	{
		_a11 += rhs._a11; _a12 += rhs._a12;
		_a21 += rhs._a21; _a22 += rhs._a22;
		return *this;
	}

	Mat22 & lshift(const size_t n)
	{
		_a11.lshift(n); _a12.lshift(n);
		_a21.lshift(n); _a22.lshift(n);
		return *this;
	}

	Mat22 & mul_right(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a12);
		_a11 *= rhs._a11; _a11.addmul(_a12, rhs._a21);
		t.addmul(_a12, rhs._a22); _a12.swap(t);

		t.mul(_a21, rhs._a12);
		_a21 *= rhs._a11; _a21.addmul(_a22, rhs._a21);
		t.addmul(_a22, rhs._a22); _a22.swap(t);

		return *this;
	}

	Mat22 & mul_right_div(const Mat22 & rhs, const gint & d)
	{
		gint t;

		t.mul(_a11, rhs._a12);
		_a11 *= rhs._a11; _a11.addmul(_a12, rhs._a21); _a11.divexact(d);
		t.addmul(_a12, rhs._a22); _a12.swap(t); _a12.divexact(d);

		t.mul(_a21, rhs._a12);
		_a21 *= rhs._a11; _a21.addmul(_a22, rhs._a21); _a21.divexact(d);
		t.addmul(_a22, rhs._a22); _a22.swap(t); _a22.divexact(d);

		return *this;
	}

	Mat22 & mul_left(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a21);
		_a11 *= rhs._a11; _a11.addmul(_a21, rhs._a12);
		t.addmul(_a21, rhs._a22); _a21.swap(t);

		t.mul(_a12, rhs._a21);
		_a12 *= rhs._a11; _a12.addmul(_a22, rhs._a12);
		t.addmul(_a22, rhs._a22); _a22.swap(t);

		return *this;
	}

	void split(Mat22 & hi, Mat22 & lo, const size_t n) const
	{
		// We may have hi.a_11/hi.a_21 > a_11/a_21 and hi.a_12/hi.a_22 < a_12/a_22 if hi.a_i = a_i >> (n * GMP_LIMB_BITS)
		// If hi.a_21 = (a_21 >> (n * GMP_LIMB_BITS)) + 1 then hi.a_11/hi.a_21 < a_11/a_21.
		// If hi.a_12 = (a_12 >> (n * GMP_LIMB_BITS)) -+1 then hi.a_12/hi.a_22 > a_12/a_22.
		// Since a_11/a_21 < alpha < a_12/a_22, we still have hi.a_11/hi.a_21 < alpha < hi.a_12/hi.a_22.
		_a11.split(hi._a11, lo._a11, n, false); _a12.split(hi._a12, lo._a12, n, true);
		_a21.split(hi._a21, lo._a21, n, true); _a22.split(hi._a22, lo._a22, n, false);
	}

	void init_gcf(const gint & N)
	{
		_a11 = 0u; _a12 = 1u;
		_a21 = N; _a21 *= 2u; _a22 = _a21;
	}

	void cf_mul(const gint & coefficient)
	{
		_a11.submul(coefficient, _a21); _a11.swap(_a21);
		_a12.submul(coefficient, _a22); _a12.swap(_a22);
	}

	void cf_mul_revert(const gint & coefficient)
	{
		_a11.swap(_a21); _a11.addmul(coefficient, _a21);
		_a12.swap(_a22); _a12.addmul(coefficient, _a22);
	}

	bool get_cf_coefficient(gint & coefficient1, gint & coefficient2)
	{
		gint t;

		coefficient1.div(_a11, _a21);
		t.div(_a12, _a22);
		if (coefficient1 != t) return false;
		cf_mul(coefficient1);

		coefficient2.div(_a11, _a21);
		t.div(_a12, _a22);
		const bool success = (coefficient2 == t);
		if (success) cf_mul(coefficient2); else cf_mul_revert(coefficient1);
		return success;
	}
};
