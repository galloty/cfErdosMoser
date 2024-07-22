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
	for (size_t j = 0; j < n / 2; ++j)
	{
		std::cout << n << ", " << j << ": " << w[j].get() << std::endl;
	}
}

static void square(Zp * const P, const Zp * const w, const size_t n)
{
	for (size_t s = 1, m = n / 2; m >= 1; s *= 2, m /= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const size_t k = 2 * j * m;
			const Zp u0 = P[k + 0 * m], u1 = P[k + 1 * m];
			P[k + 0 * m] = u0 + u1; P[k + 1 * m] = u0 - u1;

			for (size_t i = 1; i < m; ++i)
			{
				const Zp w_si = w[s * i];
				const size_t k = 2 * j * m + i;
				const Zp u0 = P[k + 0 * m], u1 = P[k + 1 * m];
				P[k + 0 * m] = u0 + u1;
				P[k + 1 * m] = (u0 - u1) * w_si;
			}
		}
	}

	for (size_t k = 0; k < n; ++k) P[k] *= P[k];

	for (size_t s = n / 2, m = 1; s >= 1; s /= 2, m *= 2)
	{
		for (size_t j = 0; j < s; ++j)
		{
			const size_t k = 2 * j * m;
			const Zp u0 = P[k + 0 * m], u1 = P[k + 1 * m];
			P[k + 0 * m] = u0 + u1; P[k + 1 * m] = u0 - u1;

			for (size_t i = 1; i < m; ++i)
			{
				const Zp w_si_inv = -w[n / 2 - s * i];
				const size_t k = 2 * j * m + i;
				const Zp u0 = P[k + 0 * m], u1 = P[k + 1 * m] * w_si_inv;
				P[k + 0 * m] = u0 + u1;
				P[k + 1 * m] = u0 - u1;
			}
		}
	}

	const Zp r = Zp::reciprocal(n);
	for (size_t k = 0; k < n; ++k) P[k] *= r;
}

int test()
{
	{
		const size_t n = 8;

		Zp w[n]; roots(w, n);
		Zp P[n]; for (size_t i = 0; i < n; ++i) P[i] = Zp(i + 1);
		square(P, w, n);

		std::cout << "Mod(1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + 6*x^5 + 7*x^6 + 8*x^7, x^8 - 1)^2" << std::endl;
		std::cout << " = 120*x^7 + 148*x^6 + 168*x^5 + 180*x^4 + 184*x^3 + 180*x^2 + 168*x + 148" << std::endl;
		for (size_t i = n - 1, j = 0; j < n; --i, ++j)
		{
			std::cout << P[i].get();
			if (i > 0) std::cout << " x^" << i << " + ";
		} 
		std::cout << std::endl;
	}

	std::cout << std::endl;

	{
		const size_t n = 16;
		Zp w[n]; roots(w, n);
		Zp P[n]; for (size_t i = 0; i < n; ++i) P[i] = Zp(i + 1);
		square(P, w, n);

		std::cout << "Mod(1 + 2*x + 3*x^2 + 4*x^3 + 5*x^4 + 6*x^5 + 7*x^6 + 8*x^7 + 9*x^8 + 10*x^9 + 11*x^10 + 12*x^11 + 13*x^12 + 14*x^13 + 15*x^14 + 16*x^15, x^16 - 1)^2" << std::endl;
		std::cout << " = 816*x^15 + 936*x^14 + 1040*x^13 + 1128*x^12 + 1200*x^11 + 1256*x^10 + 1296*x^9 + 1320*x^8 + 1328*x^7 + 1320*x^6 + 1296*x^5 + 1256*x^4 + 1200*x^3 + 1128*x^2 + 1040*x + 936" << std::endl;
		for (size_t i = n - 1, j = 0; j < n; --i, ++j)
		{
			std::cout << P[i].get();
			if (i > 0) std::cout << " x^" << i << " + ";
		} 
		std::cout << std::endl;
	}
	return EXIT_SUCCESS;
}
