/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
#include <chrono>
#include <vector>
#include <set>

#include "factor.h"
#include "mod64.h"
#include "gfloat.h"
#include "gint.h"
#include "mat22.h"
#include "pio.h"

// See Yves Gallot, Pieter Moree, Wadim Zudilin,
// The Erdős-Moser equation 1^k + 2^k + ... + (m−1)^k = m^k revisited using continued fractions,
// Math. Comp. 80 (2011), 1221-1237.

class CF
{
private:
	Heap & _heap;
	gint _N, _cond_b, _a_j;
	gfloat _q_j, _q_jm1;
	uint64_t _q_j_mod, _q_jm1_mod;
	uint64_t _j;
	Factor _factor;
	// { p : 3 is a primitive root modulo p } such that 6 * (5*7*17*19*29*31*43)^2 < 2^64
	const std::vector<uint32_t> _P_pr3 = { 5, 7, 17, 19, 29, 31, 43 };
	static const uint64_t _mod_q = 6 * uint64_t(5*7*17*19*29*31*43) * (5*7*17*19*29*31*43);

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

	void _gcf_matrix_divisor(Mat22 & M, gint & d, const uint64_t n, const uint64_t size)
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

			d = 1u;
			for (uint64_t i = 0; i < size; ++i)
			{
				const uint64_t n_i = n + i;
				const uint64_t p = _factor.smallest(n_i);
				// if n_i = p^k is a prime power then pp = p else pp = 1 (exponential of Mangoldt function).
				uint64_t m = 1; if (p != n_i) { m = n_i; while (m % p == 0) m /= p; }
				const uint64_t d_i = (m == 1) ? n_i / p : n_i;	// divides by pp
				d *= d_i;
			}
		}
		else
		{
			_gcf_matrix_divisor(M, d, n, size / 2);
			Mat22 M_r; gint d_r; _gcf_matrix_divisor(M_r, d_r, n + size / 2, size / 2);
			M.mul_right(M_r); d *= d_r;
		}
	}

	double gcf_matrix_divisor(Mat22 & M, gint & divisor, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_gcf_matrix_divisor(M, divisor, n, size);
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
		if (count1 == 0) { M.combine(M_lo, Mcf, n / 2); return 0; }

		// M_hi is Mcf * M.hi then Mcf * M = (M_hi << (n/2 * GMP_LIMB_BITS)) + Mcf * M_lo
		M_lo.mul_left(Mcf); M.combine(M_lo, Mcf, n / 2);

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
		if (count2 == 0) { M.combine(M_lo, Mcf2, n2 / 2); return count1; }

		M_lo.mul_left(Mcf2); M.combine(M_lo, Mcf2, n2 / 2);

		Mcf.mul_left(Mcf2);

		return count1 + count2;
	}

	double cf_reduce(Mat22 & M, bool & found)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		Mat22 Mcf; _cf_reduce_half(0, M, Mcf, found);

		const gfloat f11 = Mcf.get11().to_float(), f12 = Mcf.get12().to_float();
		const gfloat f21 = Mcf.get21().to_float(), f22 = Mcf.get22().to_float();
		const gfloat tf = f11 * _q_jm1 - f12 * _q_j;
		_q_j = f22 * _q_j - f21 * _q_jm1; _q_jm1 = tf;

		// Update the denominators of the regular continued fraction q_{j-1} and q_j
		const uint64_t i11 = Mcf.get11() % _mod_q, i12 = Mcf.get12() % _mod_q;
		const uint64_t i21 = Mcf.get21() % _mod_q, i22 = Mcf.get22() % _mod_q;
		Mod64 mod64(_mod_q);
		const uint64_t ti = mod64.sub(mod64.mul(i11, _q_jm1_mod), mod64.mul(i12, _q_j_mod));
		_q_j_mod = mod64.sub(mod64.mul(i22, _q_j_mod), mod64.mul(i21, _q_jm1_mod)); _q_jm1_mod = ti;

		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool condition_c() const
	{
		const uint64_t q_j = _q_jm1_mod;	// a step backward
		const uint64_t g6 = q_j % 6u;
		return ((g6 == 1) || (g6 == 5));	// (c) gcd(q_{j , 6) = 1
	}

	bool condition_d() const
	{
		// a step backward
		const uint64_t j = _j - 1;
		const gint & a_jp1 = _a_j;
		const uint64_t q_j_mod = _q_jm1_mod;
		const gfloat & q_j = _q_jm1;

		const uint64_t g6 = q_j_mod % 6u;

		std::ostringstream ss;
		ss << _N.to_string() << ", " << j << ", " << a_jp1.to_string() << ", "
		   << q_j.to_string() << ", " << ((g6 == 1) ? "+" : "-") << "1";

		bool cond_d = true;
		for (const uint32_t p : _P_pr3)
		{
			const uint64_t r = q_j_mod % (p * p);
			if ((r != 0) && (r % p == 0)) { cond_d = false; ss << ", " << p; break; }
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
		_q_j_mod = 0u; _q_jm1_mod = 1u;
		_q_j = gfloat(0, 0); _q_jm1 = gfloat(1, 0);

		const std::string Nstr = N.to_string();
		double time_gcf_matrix_divisor = 0, time_gcf_mul_div = 0, time_cf_reduce = 0, time_elapsed = 0, prev_time_elapsed = 0;
		size_t M_min_size = 0, Mgcf_size = 0;	// divisor_size = 0, M_max_size = 0;
		uint64_t j_prev = _j;

		bool found = false;
		while (!found)
		{
			{
				// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
				// Compute nstep terms of the matrix and its divisor starting at n.
				Mat22 Mgcf; gint divisor;
				time_gcf_matrix_divisor += gcf_matrix_divisor(Mgcf, divisor, n, nstep);
				Mgcf_size = Mgcf.get_byte_count();
				// divisor_size = divisor.get_byte_count();

				time_gcf_mul_div += gcf_mul_div(M, Mgcf, divisor);
			}
			n += nstep;

			// M_max_size = M.get_byte_count();

			time_cf_reduce += cf_reduce(M, found);

			M_min_size = M.get_byte_count();

			if (Mgcf_size < M_min_size / 4) nstep *= 2;

			time_elapsed = time_gcf_matrix_divisor + time_gcf_mul_div + time_cf_reduce;

			if (found) found = condition_c();

			if ((time_elapsed - prev_time_elapsed > 10) || found)
			{
				prev_time_elapsed = time_elapsed;
				std::ostringstream ss;
				ss	<< "N = " << Nstr << ", n = " << n << " [" << nstep << "], j = " << _j << " (+" << _j - j_prev << "), q_j = " << _q_j.to_string() << ", "
					// << "M_max: " << M_max_size << ", M_min: " << M_min_size << ", Mgcf: " << Mgcf_size << ", divisor: " << divisor_size << ", "
					<< "memory size: " << std::setprecision(3)
					<< _heap.get_max_size() / (1u << 20) << " + " << _heap.get_max_size_gmp() / (1u << 20) << " MB ("
					<< _heap.get_max_block_size() / (1u << 20) << " + " << _heap.get_max_block_size_gmp() / (1u << 20) << " MB), "
					<< "elapsed time: " << format_time(time_elapsed) << " ("	// "." << std::endl;
					<< time_gcf_matrix_divisor * 100 / time_elapsed << "% + "
					<< time_gcf_mul_div * 100 / time_elapsed << "% + "
					<< time_cf_reduce * 100 / time_elapsed << "%)." << std::endl;
				pio::print(ss.str());
				j_prev = _j;
			}

			if (found) found = condition_d();
		}

		_N.reset(); _cond_b.reset(); _a_j.reset();
		M.reset();
		std::cout << "Memory size: " << _heap.get_size() << " + " << _heap.get_size_gmp() << " B." << std::endl;
	}
};

int main()
{
	try
	{
		// We must have N | 2^8 * 3^5 * 5^4 * 7^3 * 11^2 * 13^2 * 17^2 * 19^2 * 23 * ... * 199,
		// where the three dots represent the product of the primes between 23 and 199. See
		// P. Moree, H. J. J. te Riele, and J. Urbanowicz,
		// Divisibility properties of integers x, k satisfying 1^ + ... + (x − 1)^k = x^k,
		// Math. Comp. 63 (1994), 799-815.

		// N_max = 2^8 * 3^5 * 5^4 * 7^3
		static const std::vector<uint32_t> fN = { 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7 };

		if (sizeof(mp_limb_t) < 8) throw std::runtime_error("32-bit computing is not supported");

		Heap heap;
		CF cf(heap);

		gint N;
		for (size_t d = 0; d < fN.size(); ++d)
		{
			N = 1u; for (size_t i = 0; i < d; ++i) N *= fN[i];
			cf.solve(N);
		}

		// Test all divisors of N_max
		std::set<uint64_t> Nset;
		for (uint64_t n2 = 1; n2 <= 256; n2 *= 2)
			for (uint64_t n3 = 1; n3 <= 243; n3 *= 3)
				for (uint64_t n5 = 1; n5 <= 625; n5 *= 5)
					for (uint64_t n7 = 1; n7 <= 343; n7 *= 7) Nset.insert(n2 * n3 * n5 * n7);
		Nset.erase(1);

		for (const uint64_t n : Nset) { N = n; cf.solve(N); }
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << "." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
