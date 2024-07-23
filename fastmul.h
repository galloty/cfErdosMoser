/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

#include <vector>
#include <iostream>

class Zp
{
private:
	static const uint64_t _p = (((uint64_t(1) << 32) - 1) << 32) + 1;	// 2^64 - 2^32 + 1
	static const uint64_t _primroot = 7;
	uint64_t _n;

public:
	Zp() {}
	explicit Zp(const uint64_t n) : _n(n) {}

	static uint64_t get_p() { return _p; }
	uint64_t get() const { return _n; }

	Zp operator-() const { return Zp((_n != 0) ? _p - _n : 0); }

	Zp & operator+=(const Zp & rhs) { const uint64_t c = (_n >= _p - rhs._n) ? _p : 0; _n += rhs._n - c; return *this; }
	Zp & operator-=(const Zp & rhs) { const uint64_t c = (_n < rhs._n) ? _p : 0; _n -= rhs._n - c; return *this; }
	Zp operator+(const Zp & rhs) const { Zp r = *this; r += rhs; return r; }
	Zp operator-(const Zp & rhs) const { Zp r = *this; r -= rhs; return r; }

	Zp & operator*=(const Zp & rhs)
	{
		const __uint128_t t = _n * __uint128_t(rhs._n);
		const uint64_t lo = uint64_t(t), hi = uint64_t(t >> 64);
		// hi.hi * 2^96 + hi.lo * 2^64 + lo = lo + hi.lo * 2^32 - (hi.hi + hi.lo)
		*this = Zp(lo) + Zp(hi << 32) - Zp((hi >> 32) + uint32_t(hi));
		return *this;
	}
	Zp operator*(const Zp & rhs) const { Zp r = *this; r *= rhs; return r; }

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

	static Zp reciprocal(const uint64_t n) { return -Zp((_p - 1) / n); }
	static const Zp primroot_n(const uint64_t n) { return Zp(_primroot).pow((_p - 1) / n); }
};

static void roots(Zp * const w, const size_t n)
{
	const Zp r = Zp::primroot_n(n);
	for (size_t j = 0; j < n / 2; ++j) w[j] = r.pow(j);
	// for (size_t j = 0; j < n / 2; ++j)
	// {
	// 	std::cout << n << ", " << j << ": " << w[j].get() << std::endl;
	// }
}

static void forward(Zp * const x, const Zp * const w, const size_t n)
{
	for (size_t s = 1, m = n / 2; m >= 1; s *= 2, m /= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const size_t k = 2 * j * m;
			const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
			x[k + 0 * m] = u0 + u1; x[k + 1 * m] = u0 - u1;

			for (size_t i = 1; i < m; ++i)
			{
				const Zp w_si = w[s * i];
				const size_t k = 2 * j * m + i;
				const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
				x[k + 0 * m] = u0 + u1;
				x[k + 1 * m] = (u0 - u1) * w_si;
			}
		}
	}
}

static void backward(Zp * const x, const Zp * const w, const size_t n)
{
	for (size_t s = n / 2, m = 1; s >= 1; s /= 2, m *= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const size_t k = 2 * j * m;
			const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
			x[k + 0 * m] = u0 + u1; x[k + 1 * m] = u0 - u1;

			for (size_t i = 1; i < m; ++i)
			{
				const Zp w_si_inv = -w[n / 2 - s * i];
				const size_t k = 2 * j * m + i;
				const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m] * w_si_inv;
				x[k + 0 * m] = u0 + u1;
				x[k + 1 * m] = u0 - u1;
			}
		}
	}

	const Zp r = Zp::reciprocal(n);
	for (size_t k = 0; k < n; ++k) x[k] *= r;
}

static void mul(Zp * const x, const Zp * const y, const size_t n)
{
	for (size_t k = 0; k < n; ++k) x[k] *= y[k];
}

static void set(Zp * const xp, const size_t n, const uint64_t * const x, const size_t x_size)
{
	for (size_t i = 0; i < x_size; ++i)
	{
		uint64_t x_i = x[i];
		for (size_t j = 0; j < 4; ++j) { xp[4 * i + j] = Zp(uint16_t(x_i)); x_i >>= 16; }
	}
	for (size_t i = 4 * x_size; i < n; ++i) xp[i] = Zp(0);
}

static void get(uint64_t * const x, const size_t x_size, const Zp * const xp)
{
	uint64_t t = 0;
	for (size_t i = 0; i < x_size; ++i)
	{
		uint64_t x_i[4];
		for (size_t j = 0; j < 4; ++j) { t += xp[4 * i + j].get(); x_i[j] = uint16_t(t); t >>= 16; }
		x[i] = x_i[0] | (x_i[1] << 16) | (x_i[2] << 32) | (x_i[3] << 48);
	}
}

// x_size >= y_size >= 16, z_size = x_size + y_size
inline void fmul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
{
	size_t e = ((x_size & (~x_size + 1)) == x_size) ? 0 : 1;	// power of two
	for (size_t r = x_size; r != 1; r /= 2) ++e;				// x_size <= 2^e

	const size_t n = size_t(1) << (e + 1 + 2);					// +1: twice the size, +2: 64-bit => 16-bit

	Zp * const w = new Zp[n];
	roots(w, n);

	Zp * const xp = new Zp[n];
	set(xp, n, x, x_size);
	Zp * const yp = new Zp[n];
	set(yp, n, y, y_size);

	forward(xp, w, n);
	forward(yp, w, n);
	mul(xp, yp, n);
	delete[] yp;
	backward(xp, w, n);
	delete[] w;

	get(z, x_size + y_size, xp);
	delete[] xp;
}
