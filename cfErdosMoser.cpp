/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <iostream>
#include <cmath>

#include <gmp.h>

class Mat22
{
private:
	mpz_t _a11, _a12;
	mpz_t _a21, _a22;

public:
	Mat22() { mpz_inits(_a11, _a12, _a21, _a22, nullptr); }
	~Mat22() { mpz_clears(_a11, _a12, _a21, _a22, nullptr); }

	void set(const mpz_t & a11, const mpz_t & a12, const mpz_t & a21, const mpz_t & a22)
	{
		mpz_set(_a11, a11); mpz_set(_a12, a12); mpz_set(_a21, a21); mpz_set(_a22, a22);
	}

	const mpz_t & get11() const { return _a11; }
	const mpz_t & get12() const { return _a12; }
	const mpz_t & get21() const { return _a21; }
	const mpz_t & get22() const { return _a22; }

	Mat22 & operator*=(const Mat22 & rhs)
	{
		mpz_t t11, t12, t21, t22; mpz_inits(t11, t12, t21, t22, nullptr);
		mpz_mul(t11, _a11, rhs._a11); mpz_addmul(t11, _a12, rhs._a21);
		mpz_mul(t12, _a11, rhs._a12); mpz_addmul(t12, _a12, rhs._a22);
		mpz_mul(t21, _a21, rhs._a11); mpz_addmul(t21, _a22, rhs._a21);
		mpz_mul(t22, _a21, rhs._a12); mpz_addmul(t22, _a22, rhs._a22);
		mpz_swap(_a11, t11); mpz_swap(_a12, t12); mpz_swap(_a21, t21); mpz_swap(_a22, t22);
		mpz_clears(t11, t12, t21, t22, nullptr);
		return *this;
	}

	void reduce()
	{
		mpz_t t; mpz_init(t);
		mpz_gcd(t, _a11, _a12); mpz_gcd(t, t, _a21); mpz_gcd(t, t, _a22);
		// mpz_out_str(stdout, 10, t); std::cout << std::endl;
		mpz_divexact(_a11, _a11, t); mpz_divexact(_a21, _a21, t); mpz_divexact(_a12, _a12, t); mpz_divexact(_a22, _a22, t);
		// mpz_gcd(t, _a11, _a21); mpz_divexact(_a11, _a11, t); mpz_divexact(_a21, _a21, t);
		// mpz_gcd(t, _a12, _a22); mpz_divexact(_a12, _a12, t); mpz_divexact(_a22, _a22, t);
		mpz_clear(t);
	}

	bool getCoefficient(mpz_t & coefficient)
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

	void print() const
	{
		mpz_out_str(stdout, 10, _a11); std::cout << std::endl;
		mpz_out_str(stdout, 10, _a21); std::cout << std::endl;
		mpz_out_str(stdout, 10, _a12); std::cout << std::endl;
		mpz_out_str(stdout, 10, _a22); std::cout << std::endl;
	}
};

int main()
{
	Mat22 M, Mg, Mn;
	mpz_t zero, one, N, q_j, q_jm1, b_n, c_n, a_j, t; mpz_inits(zero, one, N, q_j, q_jm1, b_n, c_n, a_j, t, nullptr);
	mpz_set_ui(zero, 0); mpz_set_ui(one, 1);

	for (uint64_t i = 1; i <= 256; i *= 2)
	{
		mpz_set_ui(N, i);
		mpz_set_ui(q_j, 0); mpz_set_ui(q_jm1, 1);

		mpz_add(t, N, N);
		M.set(zero, one, t, t);
		uint64_t n = 2;
		uint64_t j = 0;

		bool found = false;
		while (!found)
		{
			Mg.set(one, zero, zero, one);
			for (size_t k = 0; k < 64; ++k)
			{
				mpz_set_ui(b_n, n / 2);
				mpz_set_ui(c_n, (n % 2 == 0) ? 2 : n);
				Mn.set(zero, b_n, one, c_n);
				Mg *= Mn;
				++n;
			}

			// std::cout << mpz_sizeinbase(Mg.get11(), 2) << ", ";
			// M.print();
			Mg.reduce();
			// M.print();
			// std::cout << mpz_sizeinbase(Mg.get11(), 2) << ", ";

			M *= Mg;

			while (M.getCoefficient(a_j))
			{
				mpz_addmul(q_jm1, a_j, q_j); mpz_swap(q_j, q_jm1);
				if (j % 2 == 1)	// j - 1 is even
				{
					mpz_mul_ui(t, N, 180); mpz_sub_ui(t, t, 2);
					if (mpz_cmp(a_j, t) >= 0)	// a_j >= 180N - 2
					{
						const unsigned long g6 = mpz_mod_ui(t, q_jm1, 6);	// gcd(q_{j - 1} , 6) = 1
						if ((g6 == 1) || (g6 == 5))
						{
							long e2; const double x = mpz_get_d_2exp(&e2, q_jm1);
							const double lx = std::log10(x) + e2 * std::log10(2.0);
							const uint64_t e10 = uint64_t(lx); const double m10 = std::pow(10.0, lx - double(e10));
							std::cout << i << ", " << j - 1 << ", ";
							mpz_out_str(stdout, 10, a_j);
							std::cout << ", " << m10 << "e" << e10;
							std::cout << ", " << ((g6 == 1) ? "+" : "-") << "1";
							std::cout << std::endl;
							found = true;
							break;
						}
					}
				}
				// mpz_out_str(stdout, 10, t); std::cout << " ";
				++j;
			}
			// std::cout << " : " << j << ", ";
			// std::cout << mpz_sizeinbase(M.get11(), 2) << std::endl;
		}
	}

	// mpf_set_default_prec(1000);
	// mpf_t f, d; mpf_inits(f, d, nullptr);
	// mpf_set_z(f, M.get11()); mpf_set_z(d, M.get21()); mpf_div(f, f, d);	mpf_out_str(stdout, 10, 120, f); std::cout << std::endl;
	// mpf_set_z(f, M.get12()); mpf_set_z(d, M.get22()); mpf_div(f, f, d);	mpf_out_str(stdout, 10, 120, f); std::cout << std::endl;
	// mpf_clears(f, d, nullptr);

	mpz_clears(zero, one, N, q_j, q_jm1, b_n, c_n, a_j, t, nullptr);

	return EXIT_SUCCESS;
}
