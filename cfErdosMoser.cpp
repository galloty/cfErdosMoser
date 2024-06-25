/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <array>

#include "factor.h"
#include "gint.h"
#include "mat22.h"

class CF
{
private:
	gint _N, _cond_b, _q_j, _q_jm1, _a_j;
	std::vector<uint32_t> _PN;
	uint64_t _j;
	Factor _factor;

public:
	CF() {}
	virtual ~CF() {}

private:
	static void _gcf_matrix(Mat22 & M, const uint64_t n, const uint64_t size)
	{
		if (size == 8)
		{
			M.set_gcf(n);
			Mat22 M_i;
			for (uint64_t i = 1; i < size; ++i)
			{
				M_i.set_gcf(n + i);
				M.mul_right(M_i);
			}
		}
		else
		{
			_gcf_matrix(M, n, size / 2);
			Mat22 M_r; _gcf_matrix(M_r, n + size / 2, size / 2);
			M.mul_right(M_r);
		}
	}

	static double gcf_matrix(Mat22 & M, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_gcf_matrix(M, n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void _cfg_divisor(gint & d, const uint64_t n, const uint64_t size)
	{
		if (size == 8)
		{
			d = 1u;
			gint d_i;
			for (uint64_t i = 0; i < size; ++i)
			{
				const uint64_t n_i = n + i;
				const uint64_t p = _factor.smallest(n_i);
				// if n = p^k is a prime power then pp = p else pp = 1 (exponential of Mangoldt function).
				uint64_t m = 1; if (p != n_i) { m = n_i; while (m % p == 0) m /= p; }
				uint64_t pp = (m == 1) ? p : 1;
				d_i = n_i / pp;
				d *= d_i;
			}
		}
		else
		{
			_cfg_divisor(d, n, size / 2);
			gint d_r; _cfg_divisor(d_r, n + size / 2, size / 2);
			d *= d_r;
		}
	}

	double cfg_divisor(gint & divisor, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_cfg_divisor(divisor, n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	static double gcf_mul(Mat22 & M, const Mat22 & Mgcf)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.mul_right(Mgcf);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	static double gcf_div(Mat22 & M, const gint & divisor)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.div(divisor);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool check_conditions()
	{
		if (_j % 2 == 1)	// (a): j - 1 is even
		{
			if (_a_j >= _cond_b)	// (b) a_j >= 180N - 2
			{
				const uint32_t g6 = _q_jm1 % 6u;
				if ((g6 == 1) || (g6 == 5))		// (c) gcd(q_{j - 1} , 6) = 1
				{
					double x_2; long e_2; _q_jm1.get_d_2exp(x_2, e_2);
					const double log10_q = std::log10(x_2) + e_2 * std::log10(2.0);
					const uint64_t e10 = uint64_t(log10_q); const double m10 = std::pow(10.0, log10_q - double(e10));

					_N.out(stdout);
					std::cout << ", " << _j - 1 << ", ";
					_a_j.out(stdout);
					std::cout << ", " << m10 << "*10^" << e10;
					std::cout << ", " << ((g6 == 1) ? "+" : "-") << "1";

					bool cond_d = true;
					for (const uint32_t p : _PN)
					{
						const uint32_t nu_q_jm1 = _q_jm1.nu(p);
						if (nu_q_jm1 > 0)	// p | q_{j - 1}
						{
							// nu(3^{p−1} − 1) = 1 because p != 11, 1006003 (Mirimanoff primes)
							if (nu_q_jm1 != 1 + _N.nu(p) + 1) { cond_d = false; std::cout << ", " << p; }	// (d)
						}
					}

					std::cout << std::endl;

					if (cond_d) return true;
				}
			}
		}

		return false;
	}

	bool _cf_matrix(const Mat22 & M, Mat22 & Mcf, bool & found)
	{
		const size_t n = M.get12().get_word_count();

		if (n <= 32)
		{
			bool succeed = false;

			Mcf.set_identity();
			Mat22 Mp(M);

			// a_j is the jth coefficient of the regular continued fraction
			while (Mp.get_cf_coefficient(_a_j))
			{
				// Update the regular continued fraction
				_q_jm1.addmul(_a_j, _q_j); _q_j.swap(_q_jm1);

				Mcf.cf_mul(_a_j);

				if (check_conditions()) { found = true; break; }

				++_j;
				succeed = true;
			}

			return succeed;
		}

		Mat22 M_h; M_h.split(M, n / 2);

		if (!_cf_matrix(M_h, Mcf, found)) return false;
		if (found) return true;

		M_h = M;
		M_h.mul_left(Mcf);
		M_h.split(M_h, n / 2);

		Mat22 Mp;
		if (_cf_matrix(M_h, Mp, found))
		{
			if (!found) Mcf.mul_left(Mp);
		}

		return true;
	}

	double cf_matrix(const Mat22 & M, Mat22 & Mcf, bool & found)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_cf_matrix(M, Mcf, found);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	static double cf_mul(Mat22 & M, const Mat22 & Mcf)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.mul_left(Mcf);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

public:
	void solve(const gint & N, const std::vector<uint32_t> & PN)
	{
		_N = N; _PN = PN;
		_cond_b = N; _cond_b *= 180; _cond_b -= 2u;

		_factor.init();

		// j is the index of the convergent of the regular continued fraction.
		_j = 0;

		// n is the index of the convergent of the generalized continued fraction.
		uint64_t n = 1, nstep = 64;

		// Matrix M is the remainder of the nth convergent the generalized continued fraction of log(2) / 2N
		// after j coefficients of the regular continued fraction retrieval using the Euclidean algorithm.
		Mat22 M;

		// p_0 = 0, p_1 = 1, q_0 = q_1 = 2N.
		M.init_gcf(N);

		// Denominator of the regular continued fraction: q_j and q_{j-1}.
		_q_j = 0u; _q_jm1 = 1u;

		double time_gcf_matrix = 0, time_gcf_divisor = 0, time_gcf_mul = 0, time_gcf_div = 0, time_cf_matrix = 0, time_cf_mul = 0;
		size_t M_size = 0, Mgcf_size = 0;

		bool found = false;
		while (!found)
		{
			{
				// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
				// Compute nstep terms starting at n.
				Mat22 Mgcf; time_gcf_matrix += gcf_matrix(Mgcf, n, nstep);
				// divisor of M
				gint divisor; time_gcf_divisor += cfg_divisor(divisor, n, nstep);

				time_gcf_mul += gcf_mul(M, Mgcf);
				time_gcf_div += gcf_div(M, divisor);

				Mgcf_size = std::max(Mgcf_size, Mgcf.get_byte_count());
			}
			n += nstep;

			M_size = std::max(M_size, M.get_byte_count()); 

			Mat22 Mcf; time_cf_matrix += cf_matrix(M, Mcf, found);
			time_cf_mul += cf_mul(M, Mcf);

			if (Mgcf_size < M.get_byte_count() / 2) nstep *= 2;

			const double time_elapsed = time_gcf_matrix + time_gcf_divisor + time_gcf_mul + time_gcf_div + time_cf_matrix + time_cf_mul;
			std::cout << std::setprecision(3)
				<< "j: " << _j << ", elapsed time: " << time_elapsed << ", "
				<< "M_max: " << M_size << ", M_min: " << M.get_byte_count() << ", Mgcf: " << Mgcf_size << ", q_j: " << _q_j.get_byte_count() << ", "
				<< "gcf_matrix: " << time_gcf_matrix * 100 / time_elapsed << "%, "
				<< "gcf_divisor: " << time_gcf_divisor * 100 / time_elapsed << "%, "
				<< "gcf_mul: " << time_gcf_mul * 100 / time_elapsed << "%, "
				<< "gcf_div: " << time_gcf_div * 100 / time_elapsed << "%, "
				<< "cf_matrix: " << time_cf_matrix * 100 / time_elapsed << "%, "
				<< "cf_mul: " << time_cf_mul * 100 / time_elapsed << "%."
				<< std::endl;
		}
	}
};

int main()
{
	try
	{
		static const std::array<uint32_t, 20> fN = { 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7 };

		// P(N) where N_max = 2^8 * 3^5 * 5^4 * 7^3
		static const std::array<uint32_t, 51> PNmax = { 5, 7, 17, 19, 29, 31, 43, 101, 113, 127, 163, 197, 211, 257, 281, 379, 401, 449, 487, 631, 641, 701,
														751, 811, 1373, 1601, 2647, 2801, 3137, 4001, 4481, 7001, 13721, 16001, 17011, 18523, 22051, 28001,
														30871, 34301, 54881, 70001, 122501, 137201, 160001, 280001, 708751, 1120001, 2195201, 4167451, 5488001 };
		CF cf;
		gint N;

		// for (size_t d = 0; d < fN.size(); ++d)
		size_t d = 10;
		{
			N = 1u;
			for (size_t i = 0; i < d; ++i) N *= fN[i];

			std::vector<uint32_t> PN;	// PN is a subset of PNmax
			for (const uint32_t p : PNmax) if (N % (p - 1) == 0) PN.push_back(p);

			cf.solve(N, PN);
		}
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << "." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
