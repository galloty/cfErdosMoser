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
	static std::string format_time(const double time)
	{
		uint64_t seconds = uint64_t(time), minutes = seconds / 60, hours = minutes / 60;
		seconds -= minutes * 60; minutes -= hours * 60;

		std::stringstream ss;
		ss << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
		return ss.str();
	}

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

	static double gcf_mul_div(Mat22 & M, const Mat22 & Mgcf, const gint & divisor)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.mul_right_div(Mgcf, divisor);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	size_t _cf_reduce(Mat22 & M, Mat22 & Mcf, bool & found)
	{
		size_t count = 0;
		Mcf.set_identity();

		// a_j is the jth coefficient of the regular continued fraction
		gint a_jp1, a_jp2;
		while (M.get_cf_coefficient(a_jp1, a_jp2))
		{
			Mcf.cf_mul(a_jp1); Mcf.cf_mul(a_jp2);

			_j += 2;	// j must always be odd, now a_{j+2} is a_j
			++count;

			if (a_jp2 >= _cond_b)	// We have: (a): j-1 is even, (b) a_j >= 180N - 2
			{
				_a_j = a_jp2;
				found = true;
				break;
			}
		}

		return count;
	}

	size_t _cf_reduce_half(const size_t level, Mat22 & M, Mat22 & Mcf, bool & found)
	{
		const size_t n = M.get_min_word_count();
		if (n < 32) return _cf_reduce(M, Mcf, found);

		Mat22 M_hi, M_lo; M.split(M_hi, M_lo, n / 2);

		const size_t count1 = _cf_reduce_half(level + 1, M_hi, Mcf, found);
		if (count1 == 0) return 0;

		// M_hi is Mcf * M.hi then Mcf * M = (M_hi << (n/2 * GMP_LIMB_BITS)) + Mcf * M_lo
		M.clear(); M.swap(M_hi);
		M_lo.mul_left(Mcf); M.lshift(n / 2); M += M_lo;

		if (found) return count1;

		const size_t n2 = M.get_min_word_count();
		Mat22 Mcf2;
		if (n2 < 32)
		{
			const size_t count2 = _cf_reduce(M, Mcf2, found);
			if (count2 == 0) return count1;
			Mcf.mul_left(Mcf2);
			return count1 + count2;
		}

		M.split(M_hi, M_lo, n2 / 2);
		const size_t count2 = _cf_reduce_half(level + 1, M_hi, Mcf2, found);
		if (count2 == 0) return count1;

		M.clear(); M.swap(M_hi);
		M_lo.mul_left(Mcf2); M.lshift(n2 / 2); M += M_lo;

		Mcf.mul_left(Mcf2);

		return count1 + count2;
	}

	double cf_reduce(Mat22 & M, bool & found)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		Mat22 Mcf; _cf_reduce_half(0, M, Mcf, found);
		// Update the denominators of the regular continued fraction q_{j-1} and q_j
		gint t; t.mul(Mcf.get11(), _q_jm1); t.submul(Mcf.get12(), _q_j);
		_q_j *= Mcf.get22(); _q_j.submul(Mcf.get21(), _q_jm1);
		_q_jm1.swap(t);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool condition_c() const
	{
		const gint & q_j = _q_jm1;	// a step backward
		const uint32_t g6 = q_j % 6u;
		return ((g6 == 1) || (g6 == 5));	// (c) gcd(q_{j , 6) = 1
	}

	bool condition_d() const
	{
		// a step backward
		const uint64_t j = _j - 1;
		const gint & a_jp1 = _a_j;
		const gint & q_j = _q_jm1;

		const uint32_t g6 = q_j % 6u;

		double x_2; long e_2; q_j.get_d_2exp(x_2, e_2);
		const double log10_q = std::log10(x_2) + e_2 * std::log10(2.0);
		const uint64_t e10 = uint64_t(log10_q); const double m10 = std::pow(10.0, log10_q - double(e10));

		_N.out(stdout);
		std::cout << ", " << j << ", ";
		a_jp1.out(stdout);
		std::cout << ", " << std::setprecision(10) << m10 << "*10^" << e10;
		std::cout << ", " << ((g6 == 1) ? "+" : "-") << "1";

		bool cond_d = true;
		for (const uint32_t p : _PN)
		{
			const uint32_t nu_q_j = q_j.nu(p);
			if (nu_q_j > 0)	// p | q_j
			{
				// nu(3^{p−1} − 1) = 1 because p != 11, 1006003 (Mirimanoff primes)
				if (nu_q_j != 1 + _N.nu(p) + 1) { cond_d = false; std::cout << ", " << p; }	// (d)
			}
		}

		std::cout << std::endl;

		return cond_d;
	}

public:
	void solve(const gint & N, const std::vector<uint32_t> & PN)
	{
		_N = N; _PN = PN;
		_cond_b = N; _cond_b *= 180; _cond_b -= 2u;

		_factor.init();

		// j is the index of the convergent of the regular continued fraction.
		_j = uint64_t(-1);

		// n is the index of the convergent of the generalized continued fraction.
		uint64_t n = 1, nstep = 64;

		// Matrix M is the remainder of the nth convergent the generalized continued fraction of log(2) / 2N
		// after j coefficients of the regular continued fraction retrieval using the Euclidean algorithm.
		// j and n must be odd such that a_11/a21 < a12/a22
		Mat22 M;

		// p_0 = 0, p_1 = 1, q_0 = q_1 = 2N.
		M.init_gcf(N);

		// Denominator of the regular continued fraction: q_j and q_{j-1}.
		_q_j = 0u; _q_jm1 = 1u;

		double time_gcf_matrix = 0, time_gcf_divisor = 0, time_gcf_mul_div = 0, time_cf_reduce = 0, time_elapsed = 0, prev_time_elapsed = 0;
		size_t M_max_size = 0, M_min_size = 0, Mgcf_size = 0;	// divisor_size = 0;
		uint64_t j_prev = _j, n_prev = n;

		bool found = false;
		while (!found)
		{
			{
				// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
				// Compute nstep terms starting at n.
				Mat22 Mgcf; time_gcf_matrix += gcf_matrix(Mgcf, n, nstep);
				Mgcf_size = Mgcf.get_byte_count();
				// divisor of M
				gint divisor; time_gcf_divisor += cfg_divisor(divisor, n, nstep);
				// divisor_size = divisor.get_byte_count();

				time_gcf_mul_div += gcf_mul_div(M, Mgcf, divisor);
			}
			n += nstep;

			M_max_size = M.get_byte_count();

			time_cf_reduce += cf_reduce(M, found);

			M_min_size = M.get_byte_count();

			if (Mgcf_size < M_min_size / 2) nstep *= 2;

			time_elapsed = time_gcf_matrix + time_gcf_divisor + time_gcf_mul_div + time_cf_reduce;

			if (found) found = condition_c();

			if ((time_elapsed - prev_time_elapsed > 10) || found)
			{
				prev_time_elapsed = time_elapsed;
				std::cout << std::setprecision(3)
					<< "j = " << _j << " (+" << _j - j_prev << "), n = " << n << " (+" << n - n_prev  << "), "
					// << "M_max: " << M_max_size << ", M_min: " << M_min_size << ", Mgcf: " << Mgcf_size << ", "
					// << "divisor: " << divisor_size << ", q_j: " << _q_j.get_byte_count() << ", " << std::endl
					<< "memory size: " << M_max_size * 2 / (1u << 20) << " MB, "
					<< "elapsed time: " << format_time(time_elapsed) << ": "
					<< "gcf_matrix: " << time_gcf_matrix * 100 / time_elapsed << "%, "
					<< "gcf_divisor: " << time_gcf_divisor * 100 / time_elapsed << "%, "
					<< "gcf_mul_div: " << time_gcf_mul_div * 100 / time_elapsed << "%, "
					<< "cf_reduce: " << time_cf_reduce * 100 / time_elapsed << "%." << std::endl;
				j_prev = _j, n_prev = n;
			}

			if (found) found = condition_d();
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

		for (size_t d = 0; d < fN.size(); ++d)
		// size_t d = 10;
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
