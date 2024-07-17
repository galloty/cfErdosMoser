/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>

#include "gint.h"

class Mat22u
{
private:
	guint _a11, _a12;
	guint _a21, _a22;

public:
	Mat22u() {}
	virtual ~Mat22u() {}

	const guint & get11() const { return _a11; }
	const guint & get12() const { return _a12; }
	const guint & get21() const { return _a21; }
	const guint & get22() const { return _a22; }

	size_t get_byte_count() const { return _a11.get_byte_count() + _a12.get_byte_count() + _a21.get_byte_count() + _a22.get_byte_count(); }

	void set_gcf(const uint64_t n)
	{
		_a11 = n; _a12 = 2 * n + 1; _a12 *= n;
		_a21 = 2; _a22 = 5 * n + 2;
	}

	Mat22u & mul_right(const Mat22u & rhs)
	{
		const size_t msize = _a12.get_size() + rhs._a12.get_size();
		guint t1(msize + 1), t2(msize);

		t1.mul(_a11, rhs._a12); t2.mul(_a12, rhs._a22); t1 += t2;
		t2.mul(_a12, rhs._a21); _a12.swap(t1); t1.mul(_a11, rhs._a11); t1 += t2; _a11.swap(t1);

		t1.mul(_a21, rhs._a12); t2.mul(_a22, rhs._a22); t1 += t2;
		t2.mul(_a22, rhs._a21); _a22.swap(t1); t1.mul(_a21, rhs._a11); t1 += t2; _a21.swap(t1);

		return *this;
	}
};

class Mat22
{
private:
	gint _a11, _a12;
	gint _a21, _a22;

public:
	Mat22() {}
	virtual ~Mat22() {}

	const gint & get11() const { return _a11; }
	const gint & get12() const { return _a12; }
	const gint & get21() const { return _a21; }
	const gint & get22() const { return _a22; }

	void set_zero() { _a11 = 0; _a12 = 0; _a21 = 0; _a22 = 0; }
	void set_identity() { _a11 = 1; _a12 = 0; _a21 = 0; _a22 = 1; }

	size_t get_min_size() const { return std::min(std::min(_a11.get_size(), _a12.get_size()), std::min(_a21.get_size(), _a22.get_size())); }
	size_t get_byte_count() const { return _a11.get_byte_count() + _a12.get_byte_count() + _a21.get_byte_count() + _a22.get_byte_count(); }

	Mat22 & swap(Mat22 & rhs)
	{
		_a11.swap(rhs._a11); _a12.swap(rhs._a12);
		_a21.swap(rhs._a21); _a22.swap(rhs._a22);
		return *this;
	}

	Mat22 & div(guint & d)
	{
		int right_shift; d.div_norm(right_shift);
		guint d_inv(d.get_size() + 1); d_inv.div_invert(d);
		_a11.div_exact(d, d_inv, right_shift); _a12.div_exact(d, d_inv, right_shift);
		_a21.div_exact(d, d_inv, right_shift); _a22.div_exact(d, d_inv, right_shift);
		return *this;
	}

	Mat22 & mul_right(const Mat22u & rhs)
	{
		const size_t msize = _a12.get_size() + rhs.get12().get_size();
		gint t1(msize + 1), t2(msize);

		t1.mul(_a11, rhs.get12()); t2.mul(_a12, rhs.get22()); t1 += t2;
		t2.mul(_a12, rhs.get21()); _a12.swap(t1); t1.mul(_a11, rhs.get11()); t1 += t2; _a11.swap(t1);

		t1.mul(_a21, rhs.get12()); t2.mul(_a22, rhs.get22()); t1 += t2;
		t2.mul(_a22, rhs.get21()); _a22.swap(t1); t1.mul(_a21, rhs.get11()); t1 += t2; _a21.swap(t1);

		return *this;
	}

	Mat22 & mul_left(const Mat22 & rhs)
	{
		const size_t msize = _a12.get_size() + rhs._a12.get_size();
		gint t1(msize + 1), t2(msize);

		t1.mul(_a11, rhs._a21); t2.mul(_a21, rhs._a22); t1 += t2;
		t2.mul(_a21, rhs._a12); _a21.swap(t1); t1.mul(_a11, rhs._a11); t1 += t2; _a11.swap(t1);

		t1.mul(_a12, rhs._a21); t2.mul(_a22, rhs._a22); t1 += t2;
		t2.mul(_a22, rhs._a12); _a22.swap(t1); t1.mul(_a12, rhs._a11); t1 += t2; _a12.swap(t1);

		return *this;
	}

	void split(Mat22 & lo, const size_t n)
	{
		// We may have hi.a_11/hi.a_21 > a_11/a_21 and hi.a_12/hi.a_22 < a_12/a_22 if hi.a_i = a_i >> (n * GMP_LIMB_BITS)
		// If hi.a_21 = (a_21 >> (n * GMP_LIMB_BITS)) + 1 then hi.a_11/hi.a_21 < a_11/a_21.
		// If hi.a_12 = (a_12 >> (n * GMP_LIMB_BITS)) + 1 then hi.a_12/hi.a_22 > a_12/a_22.
		// Since a_11/a_21 < alpha < a_12/a_22, we still have hi.a_11/hi.a_21 < alpha < hi.a_12/hi.a_22.
		_a11.split(lo._a11, n); _a12.split(lo._a12, n); _a12 += 1;
		_a21.split(lo._a21, n); _a21 += 1; _a22.split(lo._a22, n);
	}

	void combine(const Mat22 & lo, const Mat22 & Ml, const size_t n)
	{
		// M -= Ml * [0 1 1 0]
		_a11 -= Ml._a12; _a11.lshift(n); _a11 += lo._a11;
		_a12 -= Ml._a11; _a12.lshift(n); _a12 += lo._a12;
		_a21 -= Ml._a22; _a21.lshift(n); _a21 += lo._a21;
		_a22 -= Ml._a21; _a22.lshift(n); _a22 += lo._a22;
	}

	void init_gcf(const guint & N)
	{
		_a11 = 0; _a12 = 1;
		_a21 = gint(N); _a21 *= 2; _a22 = _a21;
	}

	void cf_mul(const guint & coefficient)
	{
		gint t(_a22.get_size() + 2);
		t.mul(_a21, coefficient); _a11 -= t; _a11.swap(_a21);
		t.mul(_a22, coefficient); _a12 -= t; _a12.swap(_a22);
	}

	void cf_mul_revert(const guint & coefficient)
	{
		gint t(_a22.get_size() + 2);
		_a11.swap(_a21); t.mul(_a21, coefficient); _a11 += t;
		_a12.swap(_a22); t.mul(_a22, coefficient); _a12 += t;
	}

	bool get_cf_coefficient(guint & coefficient1, guint & coefficient2)
	{
		guint t;

		coefficient1.quotient(_a11.get_abs(), _a21.get_abs()); t.quotient(_a12.get_abs(), _a22.get_abs());
		if (coefficient1.cmp(t) != 0) return false;
		cf_mul(coefficient1);

		coefficient2.quotient(_a11.get_abs(), _a21.get_abs()); t.quotient(_a12.get_abs(), _a22.get_abs());
		const bool success = (coefficient2.cmp(t) == 0);
		if (success) cf_mul(coefficient2); else cf_mul_revert(coefficient1);
		return success;
	}
};
