/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

#define finline	__attribute__((always_inline))

class Mod64
{
private:
	const uint64_t _n;

public:
	Mod64(const uint64_t n) : _n(n) { }

	finline uint64_t add(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a >= _n - b) ? _n : 0;
		return a + b - c;
	}

	finline uint64_t sub(const uint64_t a, const uint64_t b) const
	{
		const uint64_t c = (a < b) ? _n : 0;
		return a - b + c;
	}

	finline uint64_t mul(const uint64_t a, const uint64_t b) const
	{
		return uint64_t((a * __uint128_t(b)) % _n);
	}
};

class Zp
{
public:
	static const uint64_t p = (((uint64_t(1) << 32) - 1) << 32) + 1;	// 2^64 - 2^32 + 1

private:
	static const uint64_t _primroot = 55;	// Such that r^{(p - 1)/4} = 2^48 and r^{(p - 1)/8} = 2^24
	uint64_t _n;

	// If a, b < p then a + b < p. If a < 2^64, b < p then a + b < 2^64.
	finline static uint64_t _add(const uint64_t a, const uint64_t b) { return a + b - ((a >= p - b) ? p : 0); }
	// If a, b < p then a - b < p. If a < 2^64, b < p then a - b < 2^64.
	finline static uint64_t _sub(const uint64_t a, const uint64_t b) { return a - b + ((a < b) ? p : 0); }
	// Output < p for any t.
	finline static uint64_t _mod(const __uint128_t t)
	{
		const uint64_t lo = uint64_t(t), hi = uint64_t(t >> 64);
		// We have 2^64 = 2^32 - 1 (mod p) and 2^96 = 2^64 - 2^32 = -1 (mod p).
		// r = hi.hi * 2^96 + hi.lo * 2^64 + lo = (hi.lo * 2^32 - hi.lo) + (lo - hi.hi)
		// hi.lo * 2^32 - hi.lo <= (2^32 - 1)^2
		// -(2^32 -1) <= lo - hi.hi <= 2^64 - 1: _sub => 0 <= lo - hi.hi <= 2^64 - 1
		// 0 <= r <= 2^64 - 1 + (2^32 - 1)^2 = 2*2^64 - 2*2^32
		// _add => 0 <= r <= 2*2^64 - 2*2^32 - (2^64 - 2^32 + 1) = 2^64 - 2^32 - 1 < p
		return _add(_sub(lo, hi >> 32), (hi << 32) - uint32_t(hi));
	}

public:
	finline Zp() {}
	finline explicit Zp(const uint64_t n) : _n(n) {}

	finline uint64_t get() const { return _n; }

	finline Zp operator+(const Zp & rhs) const { return Zp(_add(_n, rhs._n)); }
	finline Zp operator-(const Zp & rhs) const { return Zp(_sub(_n, rhs._n)); }
	finline Zp operator*(const Zp & rhs) const { return Zp(_mod(_n * __uint128_t(rhs._n))); }

	finline Zp & operator*=(const Zp & rhs) { _n = _mod(_n * __uint128_t(rhs._n)); return *this; }

	finline Zp mul_i() const { return Zp(_mod(__uint128_t(_n) << 48)); }		// i = 2^48
	finline Zp mul_sqrt_i() const { return Zp(_mod(__uint128_t(_n) << 24)); }	// sqrt(i) = 2^24
	finline Zp mul_i_sqrt_i() const { return Zp(_mod(_n * __uint128_t(((uint64_t(1) << 32) - 1) << 8))); }	// 2^72 = (2^32 - 1) * 2^8 (mod p)

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

	Zp invert() const { return pow(p - 2); }

	finline static Zp reciprocal(const int ln) { return Zp(p - ((p - 1) >> ln)); }
	static const Zp primroot_n(const int ln) { return Zp(_primroot).pow((p - 1) >> ln); }
};

typedef uint64_t v4x64 __attribute__ ((vector_size(32)));

class Zp4
{
private:
	v4x64 _v;
	static constexpr v4x64 vp = { Zp::p, Zp::p, Zp::p, Zp::p };
	static constexpr v4x64 vm32 = { uint32_t(-1), uint32_t(-1), uint32_t(-1), uint32_t(-1) };

	finline static v4x64 _add(const v4x64 & a, const v4x64 & b) { return a + b - ((a >= vp - b) & vp); }
	finline static v4x64 _sub(const v4x64 & a, const v4x64 & b) { return a - b + ((a < b) & vp); }
	finline static void _mul(const v4x64 & a, const uint64_t b, v4x64 & lo, v4x64 & hi)
	{
		__uint128_t n[4]; for (size_t i = 0; i < 4; ++i) n[i] = a[i] * __uint128_t(b);
		lo = v4x64{ uint64_t(n[0]), uint64_t(n[1]), uint64_t(n[2]), uint64_t(n[3]) };
		hi = v4x64{ uint64_t(n[0] >> 64), uint64_t(n[1] >> 64), uint64_t(n[2] >> 64), uint64_t(n[3] >> 64) };
	}
	finline static void _mul(const v4x64 & a, const v4x64 & b, v4x64 & lo, v4x64 & hi)
	{
		__uint128_t n[4]; for (size_t i = 0; i < 4; ++i) n[i] = a[i] * __uint128_t(b[i]);
		lo = v4x64{ uint64_t(n[0]), uint64_t(n[1]), uint64_t(n[2]), uint64_t(n[3]) };
		hi = v4x64{ uint64_t(n[0] >> 64), uint64_t(n[1] >> 64), uint64_t(n[2] >> 64), uint64_t(n[3] >> 64) };
	}
	finline static v4x64 _mod(const v4x64 & lo, const v4x64 & hi) { return _add(_sub(lo, hi >> 32), (hi << 32) - (hi & vm32)); }

public:
	finline Zp4() {}
	finline explicit Zp4(const v4x64 & v) : _v(v) {}
	finline explicit Zp4(const Zp & n0, const Zp & n1, const Zp & n2, const Zp & n3) : _v(v4x64{ n0.get(), n1.get(), n2.get(), n3.get() }) {}

	finline Zp4 operator+(const Zp4 & rhs) const { return Zp4(_add(_v, rhs._v)); }
	finline Zp4 operator-(const Zp4 & rhs) const { return Zp4(_sub(_v, rhs._v)); }
	finline Zp4 operator*(const Zp & rhs) const { v4x64 lo, hi; _mul(_v, rhs.get(), lo, hi); return Zp4(_mod(lo, hi)); }
	finline Zp4 operator*(const Zp4 & rhs) const { v4x64 lo, hi; _mul(_v, rhs._v, lo, hi); return Zp4(_mod(lo, hi)); }

	finline Zp4 mul_i() const { return Zp4(_mod(_v << 48, _v >> (64 - 48))); }
	finline Zp4 mul_sqrt_i() const { return Zp4(_mod(_v << 24, _v >> (64 - 24))); }
	finline Zp4 mul_i_sqrt_i() const
	{
		const v4x64 lo_40 = _v << 40, hi_40 = _v >> (64 - 40), lo_8 = _v << 8, hi_8 = _v >> (64 - 8);
		const v4x64 lo = lo_40 - lo_8, hi = hi_40 - hi_8 + (lo_40 < lo_8);
		return Zp4(_mod(lo, hi));
	}

	finline Zp4 unpacklo(const Zp4 & rhs) { return Zp4(__builtin_shuffle(_v, rhs._v, v4x64{ 0, 1, 4, 5 })); }
	finline Zp4 unpackhi(const Zp4 & rhs) { return Zp4(__builtin_shuffle(_v, rhs._v, v4x64{ 2, 3, 6, 7 }));	}

	finline static void shuffle2in(Zp4 v[4])
	{
		const Zp4 t1 = v[0].unpackhi(v[2]); v[0] = v[0].unpacklo(v[2]);
		v[2] = v[1].unpacklo(v[3]); v[3] = v[1].unpackhi(v[3]);
		v[1] = t1;
	}

	finline static void shuffle2out(Zp4 v[4])
	{
		const Zp4 t2 = v[0].unpackhi(v[1]); v[0] = v[0].unpacklo(v[1]);
		v[1] = v[2].unpacklo(v[3]); v[3] = v[2].unpackhi(v[3]);
		v[2] = t2;
	}
};
