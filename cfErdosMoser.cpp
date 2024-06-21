/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <cmath>

#include <gmp.h>

#include "arith64.h"

class Mat22
{
private:
	mpz_t _a11, _a12;
	mpz_t _a21, _a22;

public:
	Mat22() { mpz_inits(_a11, _a12, _a21, _a22, nullptr); }
	~Mat22() { mpz_clears(_a11, _a12, _a21, _a22, nullptr); }

	Mat22 & operator*=(const Mat22 & rhs)
	{
		mpz_t t; mpz_init(t);

		mpz_mul(t, _a11, rhs._a12);
		mpz_mul(_a11, _a11, rhs._a11); mpz_addmul(_a11, _a12, rhs._a21);
		mpz_addmul(t, _a12, rhs._a22); mpz_swap(_a12, t);

		mpz_mul(t, _a21, rhs._a12);
		mpz_mul(_a21, _a21, rhs._a11); mpz_addmul(_a21, _a22, rhs._a21);
		mpz_addmul(t, _a22, rhs._a22); mpz_swap(_a22, t);

		mpz_clear(t);

		return *this;
	}

	Mat22 & operator/=(const mpz_t & d)
	{
		if (mpz_divisible_p(_a11, d)) mpz_divexact(_a11, _a11, d); else { std::cout << "Error" << std::endl; throw; }
		if (mpz_divisible_p(_a12, d)) mpz_divexact(_a12, _a12, d); else { std::cout << "Error" << std::endl; throw; }
		if (mpz_divisible_p(_a21, d)) mpz_divexact(_a21, _a21, d); else { std::cout << "Error" << std::endl; throw; }
		if (mpz_divisible_p(_a22, d)) mpz_divexact(_a22, _a22, d); else { std::cout << "Error" << std::endl; throw; }
		return *this;
	}

	void init_gcf(const mpz_t & N)
	{
		mpz_set_ui(_a11, 0); mpz_set_ui(_a12, 1);
		mpz_set(_a21, N); mpz_mul_2exp(_a21, _a21, 1); mpz_set(_a22, _a21);
	}

	void eval_gcf(const uint64_t n, const uint64_t size)
	{
		if (size == 1)
		{
			mpz_set_ui(_a11, 0); mpz_set_ui(_a12, n / 2);
			mpz_set_ui(_a21, 1); mpz_set_ui(_a22, (n % 2 == 0) ? 2 : n);
		}
		else
		{
			eval_gcf(n, size / 2);
			Mat22 M; M.eval_gcf(n + size / 2, size / 2);
			*this *= M;
		}
	}

	bool get_cf_coefficient(mpz_t & coefficient)
	{
		mpz_t t; mpz_init(t);
		mpz_fdiv_q(coefficient, _a11, _a21);
		mpz_fdiv_q(t, _a12, _a22);
		const bool success = (mpz_cmp(coefficient, t) == 0);
		mpz_clear(t);
		if (success)
		{
			mpz_submul(_a11, coefficient, _a21); mpz_swap(_a11, _a21);
			mpz_submul(_a12, coefficient, _a22); mpz_swap(_a12, _a22);
		}
		return success;
	}
};

static void eval_mdivisor(const uint64_t n, const uint64_t size, mpz_t & mdivisor)
{
	if (size == 1)
	{
		mpz_set_ui(mdivisor, n / isprimepower(n));
	}
	else
	{
		eval_mdivisor(n, size / 2, mdivisor);
		mpz_t rdivisor; mpz_init(rdivisor); eval_mdivisor(n + size / 2, size / 2, rdivisor);
		mpz_mul(mdivisor, mdivisor, rdivisor);
		mpz_clear(rdivisor);
	}
}

int main()
{
	static const uint32_t fN[20] = { 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7 };

	for (size_t d = 0; d < 10; ++d)
	{
		mpz_t N; mpz_init_set_ui(N, 1);
		for (size_t i = 0; i < d; ++i) mpz_mul_ui(N, N, fN[i]);
		mpz_t cond_b; mpz_init(cond_b); mpz_mul_ui(cond_b, N, 180); mpz_sub_ui(cond_b, cond_b, 2);

		// n is the index of the convergent of the generalized continued fraction,
		// j is the index of the convergent of the regular continued fraction.
		uint64_t n = 2, j = 0;

		// Matrix M is the remainder of the nth convergent the generalized continued fraction of log(2) / 2N
		// after j coefficients of the regular continued fraction retrieval using the Euclidean algorithm.
		Mat22 M;

		// p_0 = 0, p_1 = 1, q_0 = q_1 = 2N.
		M.init_gcf(N);

		// Denominator of the regular continued fraction: q_j and q_{j-1}.
		mpz_t q_j, q_jm1; mpz_init_set_ui(q_j, 0); mpz_init_set_ui(q_jm1, 1);

		mpz_t mdivisor, a_j; mpz_inits(mdivisor, a_j, nullptr);

		bool found = false;
		while (!found)
		{
			// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
			// Compute 64 terms starting at n.
			Mat22 Mg; Mg.eval_gcf(n, 64);
			eval_mdivisor(n / 2, 64 / 2, mdivisor);
			n += 64;

			M *= Mg;
			M /= mdivisor;

			// a_j is the jth coefficient of the regular continued fraction
			while (M.get_cf_coefficient(a_j))
			{
				// Update the regular continued fraction
				mpz_addmul(q_jm1, a_j, q_j); mpz_swap(q_j, q_jm1);

				if (j % 2 == 1)	// j - 1 is even
				{
					if (mpz_cmp(a_j, cond_b) >= 0)	// a_j >= 180N - 2
					{
						const unsigned long g6 = mpz_fdiv_ui(q_jm1, 6);	// gcd(q_{j - 1} , 6) = 1
						if ((g6 == 1) || (g6 == 5))
						{
							long e2; const double x = mpz_get_d_2exp(&e2, q_jm1);
							const double lx = std::log10(x) + e2 * std::log10(2.0);
							const uint64_t e10 = uint64_t(lx); const double m10 = std::pow(10.0, lx - double(e10));

							mpz_out_str(stdout, 10, N);
							std::cout << ", " << j - 1 << ", ";
							mpz_out_str(stdout, 10, a_j);
							std::cout << ", " << m10 << "*10^" << e10;
							std::cout << ", " << ((g6 == 1) ? "+" : "-") << "1";
							std::cout << std::endl;
							found = true;
							break;
						}
					}
				}

				++j;
			}
		}

		mpz_clears(N, cond_b, q_j, q_jm1, mdivisor, a_j, nullptr);
	}

	return EXIT_SUCCESS;
}
