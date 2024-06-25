/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

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

	void set_identity() { _a11 = 1u; _a12 = 0u; _a21 = 0u; _a22 = 1u; }

	void set_gcf(const uint64_t n)
	{
		_a11 = n; _a12 = 2 * n + 1; _a12 *= _a11;
		_a21 = 2u; _a22 = 5 * n + 2;
	}

	size_t get_byte_count() const { return _a11.get_byte_count() + _a12.get_byte_count() + _a21.get_byte_count() + _a22.get_byte_count(); }

	void mul_right(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a12);
		_a11 *= rhs._a11; _a11.addmul(_a12, rhs._a21);
		t.addmul(_a12, rhs._a22); _a12.swap(t);

		t.mul(_a21, rhs._a12);
		_a21 *= rhs._a11; _a21.addmul(_a22, rhs._a21);
		t.addmul(_a22, rhs._a22); _a22.swap(t);
	}

	void mul_left(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a21);
		_a11 *= rhs._a11; _a11.addmul(_a21, rhs._a12);
		t.addmul(_a21, rhs._a22); _a21.swap(t);

		t.mul(_a12, rhs._a21);
		_a12 *= rhs._a11; _a12.addmul(_a22, rhs._a12);
		t.addmul(_a22, rhs._a22); _a22.swap(t);
	}

	void div(const gint & d) { _a11.divexact(d); _a12.divexact(d); _a21.divexact(d); _a22.divexact(d); }

	void split(const Mat22 & rhs, const size_t n)
	{
		_a11.rshift(rhs._a11, n);
		_a12.rshift(rhs._a12, n); _a12 += 1;
		_a21.rshift(rhs._a21, n); _a21 += 1;
		_a22.rshift(rhs._a22, n);
	}

	void init_gcf(const gint & N)
	{
		_a11 = 0u; _a12 = 1u;
		_a21 = N; _a21 += _a21; _a22 = _a21;
	}

	void cf_mul(const gint & coefficient)
	{
		_a11.submul(coefficient, _a21); _a11.swap(_a21);
		_a12.submul(coefficient, _a22); _a12.swap(_a22);
	}

	bool get_cf_coefficient(gint & coefficient)
	{
		coefficient.div(_a11, _a21);
		gint t; t.div(_a12, _a22);
		const bool success = (coefficient == t);
		if (success) cf_mul(coefficient);
		return success;
	}
};
