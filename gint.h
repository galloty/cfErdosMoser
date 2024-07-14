/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

#include "guint.h"
#include "gfloat.h"

// giant integer
class gint
{
private:
	guint _u;
	bool _is_positive;

private:
	int _cmp(const gint & rhs) const
	{
		const bool is_positive = _is_positive, ris_positive = rhs._is_positive;
		const bool is_zero = _u.is_zero(), ris_zero = rhs._u.is_zero();

		// 0 may be positive or negative
		if (is_zero)
		{
			if (ris_zero) return 0;
			return ris_positive ? -1 : 1;
		}
		if (ris_zero) return is_positive ? 1 : -1;

		if (is_positive != ris_positive) return is_positive ? 1 : -1;
		const int ucmp = _u.cmp(rhs._u);
		return is_positive ? ucmp : -ucmp;
	}

public:
	gint() : _is_positive(true) {}
	gint(const uint64_t n) : _u(n), _is_positive(true) {}
	gint(const gint & rhs) : _u(rhs._u), _is_positive(rhs._is_positive) {}
 	virtual ~gint() {}

 	size_t get_word_count() const { return _u.get_word_count(); }
 	size_t get_byte_count() const { return _u.get_byte_count(); }

	gint & operator=(const uint64_t n) { _u = n; _is_positive = true; return *this; }
	gint & operator=(const gint & rhs) { if (&rhs != this) { _u = rhs._u; _is_positive = rhs._is_positive; } return *this; }

	gint & swap(gint & rhs) { _u.swap(rhs._u); std::swap(_is_positive, rhs._is_positive); return *this; }

 	bool operator==(const gint & rhs) const { return (_cmp(rhs) == 0); }
 	bool operator!=(const gint & rhs) const { return (_cmp(rhs) != 0); }
 	bool operator>=(const gint & rhs) const { return (_cmp(rhs) >= 0); }

	gint & operator+=(const uint64_t n) { if (_is_positive) _u += n; else _is_positive = !_u.sub(n); return *this; }
	gint & operator-=(const uint64_t n) { if (_is_positive) _is_positive = _u.sub(n); else _u += n; return *this; }
	gint & operator*=(const uint64_t n) { _u *= n; return *this; }

 	static void mod_init(const uint64_t n, uint64_t & n_inv, int & n_e) { guint::mod_init(n, n_inv, n_e); }

	uint64_t mod(const uint64_t n, const uint64_t n_inv, const int n_e, const uint64_t n_f) const
	{
		const uint64_t remainder = _u.mod(n, n_inv, n_e, n_f);
		if (remainder == 0) return 0;
		return _is_positive ? remainder : n - remainder;
	}

	gint & operator+=(const gint & rhs)
	{
		if (_is_positive == rhs._is_positive) _u += rhs._u;
		else _is_positive = _is_positive ? _u.sub(rhs._u) : !_u.sub(rhs._u);
		return *this;
	}

	gint & operator-=(const gint & rhs)
	{
		if (_is_positive != rhs._is_positive) _u += rhs._u;
		else _is_positive = _is_positive ? _u.sub(rhs._u) : !_u.sub(rhs._u);
		return *this;
	}

	gint & mul(const gint & x, const gint & y)	// *this != x, *this != y
	{
		_u.mul(x._u, y._u);
		_is_positive = (x._is_positive == y._is_positive);
		return *this;
	}

	gint & operator*=(const gint & rhs) { gint t; t.mul(*this, rhs); swap(t); return *this; }

 	gint & lshift(const size_t n) { _u.lshift(n); return *this; }
 	gint & rshift(const size_t n) { _u.rshift(n); return *this; }

	gint & div_norm(int & right_shift) { _u.div_norm(right_shift); return *this; }

	// 2^{2*d_bit_size} / d: |d_inv| > |d|
	gint & div_invert(const gint & d)	// *this != d
	{
		_u.div_invert(d._u);
		_is_positive = d._is_positive;
		return *this;
	}

	gint & div_exact(const gint & d, const gint & d_inv, const int right_shift)
	{
		_u.div_exact(d._u, d_inv._u, right_shift);
		return *this;
	}

	gint & quotient(const gint & x, const gint & y)	// *this != x, *this != y
	{
		if (!x._is_positive || !y._is_positive) throw std::runtime_error("divide: negative input");
		_u.quotient(x._u, y._u);
		return *this;
	}

	void split(gint & lo, const size_t n)
	{
		_u.split(lo._u, n);
		lo._is_positive = _is_positive;
	}

	gfloat to_float() const
	{
		gfloat uf = _u.to_float();
		return _is_positive ? uf : -uf;
	}

	std::string to_string() const
	{
		const std::string str_pos = _u.to_string();
		return _is_positive ? str_pos : (std::string("-") + str_pos);
	}
};
