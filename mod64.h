/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

class Mod64
{
private:
	const uint64_t _n;

public:
	Mod64(const uint64_t n) : _n(n) { }

	uint64_t add(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a >= _n - b) ? _n : 0;
		return a + b - c;
	}

	uint64_t sub(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a < b) ? _n : 0;
		return a - b + c;
	}

	uint64_t mul(const uint64_t a, const uint64_t b) const
	{
		return uint64_t((a * __uint128_t(b)) % _n);
	}
};

class Zp
{
private:
	static const uint64_t _p = (((uint64_t(1) << 32) - 1) << 32) + 1;	// 2^64 - 2^32 + 1
	static const uint64_t _mp = (uint64_t(1) << 32) - 1;				// -p = 2^32 - 1
	static const uint64_t _primroot = 7;
	uint64_t _n;

	Zp & _mod_p(const __uint128_t t)
	{
		const uint64_t lo = uint64_t(t), hi = uint64_t(t >> 64);
		// hi.hi * 2^96 + hi.lo * 2^64 + lo = lo + hi.lo * 2^32 - (hi.hi + hi.lo)
		*this = Zp().set(lo) + Zp(hi << 32) - Zp((hi >> 32) + uint32_t(hi));
		return *this;
	}

public:
	Zp() {}
	explicit Zp(const uint64_t n) : _n(n) {}

	uint64_t get() const { return _n; }
	Zp & set(const uint64_t n) { const uint64_t c = (n >= _p) ? _mp : 0; _n = n + c; return *this; }

	Zp operator-() const { return Zp((_n != 0) ? _p - _n : 0); }

	Zp & operator+=(const Zp & rhs) { const uint64_t c = (_n >= _p - rhs._n) ? _mp : 0; _n += rhs._n + c; return *this; }
	Zp & operator-=(const Zp & rhs) { const uint64_t c = (_n < rhs._n) ? _mp : 0; _n -= rhs._n + c; return *this; }
	Zp operator+(const Zp & rhs) const { Zp r = *this; r += rhs; return r; }
	Zp operator-(const Zp & rhs) const { Zp r = *this; r -= rhs; return r; }

	Zp & operator*=(const Zp & rhs) { _mod_p(_n * __uint128_t(rhs._n)); return *this; }
	Zp operator*(const Zp & rhs) const { Zp r = *this; r *= rhs; return r; }

	Zp mul_i() const { return Zp()._mod_p(__uint128_t(_n) << 48); }	// i = 2^48

	Zp pow(const uint64_t e) const
	{
		if (e == 0) return Zp(1);

		Zp r = Zp(1), y = *this;
		for (uint64_t i = e; i != 1; i /= 2)
		{
			if (i % 2 != 0) r *= y;
			y *= y;
		}
		r *= y;

		return r;
	}

	Zp invert() const { return pow(_p - 2); }

	static Zp reciprocal(const uint64_t n) { return -Zp((_p - 1) / n); }
	static const Zp primroot_n(const uint64_t n) { return Zp(_primroot).pow((_p - 1) / n); }
};
