/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>

inline uint64_t mulmod(const uint64_t x, const uint64_t y, const uint64_t m)
{
	return uint64_t((x * __uint128_t(y)) % m);
}

inline uint64_t powmod(const uint64_t x, const uint64_t n, const uint64_t m)
{
	if (n == 0) return 1;

	uint64_t r = 1, y = x;
	for (uint64_t i = n; i != 1; i /= 2)
	{
		if (i % 2 != 0) r = mulmod(r, y, m);
		y = mulmod(y, y, m);
	}
	return mulmod(r, y, m);
}

inline uint64_t gcd(const uint64_t x, const uint64_t y)
{
	uint64_t a = x, b = y;
	while (b != 0)
	{
		const uint64_t t = b;
		b = a % b;
		a = t;
	}
	return a;
}

inline bool spsp(const uint64_t n, const uint64_t p)
{
	// n - 1 = 2^k * r
	uint64_t r = n - 1;
	int k = 0;
	for (; r % 2 == 0; r /= 2) ++k;

	uint64_t x = powmod(p, r, n);
	if (x == 1) return true;

	// Compute x^(2^i) for 0 <= i < n.  If any are -1, n is a p-spsp.
	for (; k > 0; --k)
	{
		if (x == n - 1) return true;
		x = mulmod(x, x, n);
	}

	return false;
}

inline bool isprime(const uint64_t n)
{
	if (n < 2) return false;
	if (n % 2 == 0) return (n == 2);
	if (n < 9) return true;

	// see https://oeis.org/A014233

	if (!spsp(n, 2)) return false;
	if (n < 2047ull) return true;

	if (!spsp(n, 3)) return false;
	if (n < 1373653ull) return true;

	if (!spsp(n, 5)) return false;
	if (n < 25326001ull) return true;

	if (!spsp(n, 7)) return false;
	if (n < 3215031751ull) return true;

	if (!spsp(n, 11)) return false;
	if (n < 2152302898747ull) return true;

	if (!spsp(n, 13)) return false;
	if (n < 3474749660383ull) return true;

	if (!spsp(n, 17)) return false;
	if (n < 341550071728321ull) return true;

	if (!spsp(n, 19)) return false;
	// if (n < 341550071728321ull) return true;

	if (!spsp(n, 23)) return false;
	if (n < 3825123056546413051ull) return true;

	if (!spsp(n, 29)) return false;
	// if (n < 3825123056546413051ull) return true;

	if (!spsp(n, 31)) return false;
	// if (n < 3825123056546413051ull) return true;

	if (!spsp(n, 37)) return false;
	return true;	// 318665857834031151167461
}

// if n = p^k is a prime power then return p else 1 (exponential of Mangoldt function).
inline uint64_t isprimepower(const uint64_t n)
{
	if (n < 4) return n;
	if ((n & (~n + 1)) == n) return 2;	// power of two
	if (n % 2 == 0) return 1;

	// see Henri Cohen, A Course in Computational Algebraic Number Theory.
	// Graduate Texts in Mathematics. Vol. 138. Springer-Verlag. Algorithm 1.7.4

	for (uint64_t a = 2; true; ++a)
	{
		const uint64_t b = powmod(a, n, n);
		const uint64_t p = gcd(n, (b >= a) ? b - a : b - a + n);
		if (p == 1) return 1;
		if (isprime(p))
		{
			uint64_t m = n; while (m % p == 0) m /= p;
			return (m == 1) ? p : 1;
		}
	}
}