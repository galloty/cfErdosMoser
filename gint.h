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
#include "checkpoint.h"

// giant integer
class gint
{
private:
	guint _u;
	bool _is_positive;

public:
	gint(const size_t alloc_size = 1) : _u(alloc_size), _is_positive(true) {}
	gint(const gint & rhs) : _u(rhs._u), _is_positive(rhs._is_positive) {}
	explicit gint(const guint & rhs) : _u(rhs), _is_positive(true) {}
 	virtual ~gint() {}

 	size_t get_size() const { return _u.get_size(); }
 	size_t get_byte_count() const { return _u.get_byte_count(); }
	bool is_positive() const { return (_is_positive || (get_size() == 0)); }
	bool is_negative() const { return (!_is_positive || (get_size() == 0)); }

	const guint & get_abs() const { return _u; }

	gint & operator=(const uint64_t n) { _u = n; _is_positive = true; return *this; }
	gint & operator=(const gint & rhs) { if (&rhs != this) { _u = rhs._u; _is_positive = rhs._is_positive; } return *this; }

	void clear() { _u.clear(); }

	gint & swap(gint & rhs) { _u.swap(rhs._u); std::swap(_is_positive, rhs._is_positive); return *this; }

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

	gint & mul(const gint & x, const guint & y)	// *this != x
	{
		_u.mul(x._u, y);
		_is_positive = x._is_positive;
		return *this;
	}

	gint & lshift(const size_t n) { _u.lshift(n); return *this; }
	gint & rshift(const size_t n) { _u.rshift(n); return *this; }

	gint & div_exact(const guint & d, const guint & d_inv, const int right_shift) { _u.div_exact(d, d_inv, right_shift); return *this; }

	void split(gint & lo, const size_t n) { _u.split(lo._u, n); lo._is_positive = _is_positive; }

	gfloat to_float() const
	{
		const gfloat uf = _u.to_float();
		return _is_positive ? uf : -uf;
	}

	std::string to_string() const
	{
		const std::string str_pos = _u.to_string();
		return _is_positive ? str_pos : (std::string("-") + str_pos);
	}

	bool read(Checkpoint & checkpoint)
	{
		uint32_t is_positive; if (!checkpoint.read(reinterpret_cast<char *>(&is_positive), sizeof(is_positive))) return false;
		_is_positive = (is_positive != 0);
		return _u.read(checkpoint);
	}

	bool write(Checkpoint & checkpoint) const
	{
		const uint32_t is_positive = _is_positive ? 1 : 0;
		if (!checkpoint.write(reinterpret_cast<const char *>(&is_positive), sizeof(is_positive))) return false;
		return _u.write(checkpoint);
	}
};
