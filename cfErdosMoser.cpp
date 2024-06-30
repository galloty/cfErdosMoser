/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <vector>
#include <set>

#include "factor.h"
#include "gint.h"
#include "mat22.h"
#include "pio.h"

// See Yves Gallot, Pieter Moree, Wadim Zudilin,
// The Erdős-Moser equation 1^k + 2^k + ... + (m−1)^k = m^k revisited using continued fractions,
// Math. Comp. 80 (2011), 1221-1237.


class CF
{
private:
	// We must have N | 2^8 * 3^5 * 5^4 * 7^3 * 11^2 * 13^2 * 17^2 * 19^2 * 23 * ... * 199,
	// where the three dots represent the product of the primes between 23 and 199. See
	// P. Moree, H. J. J. te Riele, and J. Urbanowicz,
	// Divisibility properties of integers x, k satisfying 1^ + ... + (x − 1)^k = x^k,
	// Math. Comp. 63 (1994), 799-815.

	// Here N is a divisor of N_max = 2^8 * 3^5 * 5^4 * 7^3

private:
	// { p : p < 1000 and p - 1 | N_max and 3 is not a primitive root modulo p }
	const std::vector<uint32_t> _P1_Nmax = { 11, 13, 37, 41, 61, 71, 73, 97, 109, 151, 181, 193, 241, 251, 271, 337,
											421, 433, 491, 541, 577, 601, 673, 757, 769, 883 };

	// { p : p < 1000 and 3 is a primitive root modulo p }
	const std::vector<uint32_t> _P_pr3 = { 5, 7, 17, 19, 29, 31, 43, 53, 79, 89, 101, 113, 127, 137, 139, 149, 163, 173, 197,
											199, 211, 223, 233, 257, 269, 281, 283, 293, 317, 331, 353, 379, 389, 401, 449, 461,
											463, 487, 509, 521, 557, 569, 571, 593, 607, 617, 631, 641, 653, 677, 691, 701, 739,
											751, 773, 797, 809, 811, 821, 823, 857, 859, 881, 907, 929, 941, 953, 977 };
private:
	Heap & _heap;
	gint _N, _cond_b, _q_j, _q_jm1, _a_j;
	std::vector<uint32_t> _P1_N;
	uint64_t _j;
	Factor _factor;

public:
	CF(Heap & heap) : _heap(heap) {}
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
		if (size == 32)
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
		if (size == 32)
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

		Mat22 M_lo; M.split(M_lo, n / 2);

		const size_t count1 = _cf_reduce_half(level + 1, M, Mcf, found);
		if (count1 == 0) { M.lshift(n / 2); M += M_lo; return 0; }

		// M_hi is Mcf * M.hi then Mcf * M = (M_hi << (n/2 * GMP_LIMB_BITS)) + Mcf * M_lo
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

		M.split(M_lo, n2 / 2);
		const size_t count2 = _cf_reduce_half(level + 1, M, Mcf2, found);
		if (count2 == 0) { M.lshift(n2 / 2); M += M_lo; return count1; }

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

		std::ostringstream ss;
		ss << _N.to_string() << ", " << j << ", " << a_jp1.to_string();
		ss << ", " << std::setprecision(10) << m10 << "*10^" << e10;
		ss << ", " << ((g6 == 1) ? "+" : "-") << "1";

		bool cond_d = true;
		for (const uint32_t p : _P1_N)
		{
			const uint32_t r = q_j % (p * p);
			if ((r != 0) && (r % p == 0)) { cond_d = false; ss << ", " << p; break; }
		}
		if (cond_d)
		{
			for (const uint32_t p : _P_pr3)
			{
				const uint32_t r = q_j % (p * p);
				if ((r != 0) && (r % p == 0)) { cond_d = false; ss << ", " << p; break; }
			}
		}
		ss << std::endl;
		pio::print(ss.str(), true);

		return cond_d;
	}

public:
	void solve(const gint & N)
	{
		_heap.reset_max_size();

		_N = N;

		// P1_N is a subset of P1_Nmax
		for (const uint32_t p : _P1_Nmax) if (N % (p - 1) == 0) _P1_N.push_back(p);

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
		size_t M_min_size = 0, Mgcf_size = 0;	// divisor_size = 0, M_max_size = 0;
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

			// M_max_size = M.get_byte_count();

			time_cf_reduce += cf_reduce(M, found);

			M_min_size = M.get_byte_count();

			if (Mgcf_size < M_min_size / 2) nstep *= 2;

			time_elapsed = time_gcf_matrix + time_gcf_divisor + time_gcf_mul_div + time_cf_reduce;

			if (found) found = condition_c();

			if ((time_elapsed - prev_time_elapsed > 10) || found)
			{
				prev_time_elapsed = time_elapsed;
				std::ostringstream ss;
				ss << std::setprecision(3)
					<< "j = " << _j << " (+" << _j - j_prev << "), n = " << n << " (+" << n - n_prev  << "), "
					// << "M_max: " << M_max_size << ", M_min: " << M_min_size << ", Mgcf: " << Mgcf_size << ", "
					// << "divisor: " << divisor_size << ", q_j: " << _q_j.get_byte_count() << ", " << std::endl
					<< "memory size: " << _heap.get_max_size() * 2 / (1u << 20) << " MB, "
					<< "elapsed time: " << format_time(time_elapsed) << ": "
					<< "gcf_matrix: " << time_gcf_matrix * 100 / time_elapsed << "%, "
					<< "gcf_divisor: " << time_gcf_divisor * 100 / time_elapsed << "%, "
					<< "gcf_mul_div: " << time_gcf_mul_div * 100 / time_elapsed << "%, "
					<< "cf_reduce: " << time_cf_reduce * 100 / time_elapsed << "%." << std::endl;
				pio::print(ss.str());
				j_prev = _j, n_prev = n;
			}

			if (found) found = condition_d();
		}

		_N.reset(); _cond_b.reset(); _q_j.reset(); _q_jm1.reset(); _a_j.reset();
	}
};

int main()
{
	try
	{
		// N_max = 2^8 * 3^5 * 5^4 * 7^3
		static const std::vector<uint32_t> fN = { 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7 };

		Heap heap;
		CF cf(heap);

		std::set<uint64_t> Nset;
		for (uint64_t n2 = 1; n2 <= 256; n2 *= 2)
			for (uint64_t n3 = 1; n3 <= 243; n3 *= 3)
				for (uint64_t n5 = 1; n5 <= 625; n5 *= 5)
					for (uint64_t n7 = 1; n7 <= 343; n7 *= 7) Nset.insert(n2 * n3 * n5 * n7);
		Nset.erase(1);

		gint N;
		for (const uint64_t n : Nset) { N = n; cf.solve(N); }
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << "." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
