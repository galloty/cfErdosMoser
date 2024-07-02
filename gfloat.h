/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

// giant float
class gfloat
{
private:
	// mantissa * 2^exponent
	double _mantissa;
	size_t _exponent;

	void _norm()
	{
		if (_mantissa == 0) { _exponent = 0; return; }
		if (_mantissa > 0)
		{
			while (_mantissa < 0.5) { _mantissa *= 2; --_exponent; if (_exponent == 0) break; }
			while (_mantissa >= 1) { _mantissa *= 0.5; ++_exponent; }
		}
		else
		{
			while (_mantissa > -0.5) { _mantissa *= 2; --_exponent; if (_exponent == 0) break; }
			while (_mantissa <= -1) { _mantissa *= 0.5; ++_exponent; }
		}
	}

public:
	gfloat() {}
	gfloat(const double mantissa, const size_t exponent) : _mantissa(mantissa), _exponent(exponent) {}
	virtual ~gfloat() {}

	double get_mantissa() const { return _mantissa; }
	size_t get_exponent() const { return _exponent; }

	gfloat & operator+=(const gfloat & rhs)
	{
		if (_exponent > rhs._exponent + 64) return *this;
		if (rhs._exponent > _exponent + 64) { *this = rhs; return *this; }

		int e = 0;
		if (_exponent > rhs._exponent) e = int(rhs._exponent - _exponent);
		if (_exponent < rhs._exponent) e = -int(_exponent - rhs._exponent);
		_mantissa += std::ldexpl(rhs._mantissa, e);
		_norm();
		return *this;
	}

	gfloat & operator-=(const gfloat & rhs) { *this += gfloat(-rhs._mantissa, rhs._exponent); return *this; }
	gfloat & operator*=(const gfloat & rhs) { _mantissa *= rhs._mantissa; _exponent += rhs._exponent; _norm(); return *this; }

	gfloat operator+(const gfloat & rhs) const { gfloat r = *this; r += rhs; return r; }
	gfloat operator-(const gfloat & rhs) const { gfloat r = *this; r -= rhs; return r; }
	gfloat operator*(const gfloat & rhs) const { gfloat r = *this; r *= rhs; return r; }

	// convert to base 10
	gfloat to_base10() const
	{
		if (_mantissa == 0) return gfloat(0, 0);
		const int sgn = (_mantissa < 0) ? -1 : 1;
		const double m = std::fabs(_mantissa);
		const long double lg10 = std::log10l(m) + _exponent * std::log10l(2.0L);
		size_t exponent10 = uint64_t(lg10); double mantissa10 = std::pow(10.0, double(lg10 - exponent10));
		while (mantissa10 >= 10) { mantissa10 *= 0.1; ++exponent10; }
		while (mantissa10 < 1) { mantissa10 *= 10; --exponent10; }
		return gfloat(sgn * mantissa10, exponent10);
	}

	std::string to_string() const
	{
		std::ostringstream ss;
		ss << std::setprecision(10) << _mantissa << "e+" << _exponent;
		return ss.str();
	}
};
