/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <sys/stat.h>

#include "factor.h"
#include "mod64.h"
#include "checkpoint.h"
#include "gfloat.h"
#include "heap.h"
#include "guint.h"
#include "gint.h"
#include "mat22.h"

// See Yves Gallot, Pieter Moree, Wadim Zudilin,
// The Erdős-Moser equation 1^k + 2^k + ... + (m−1)^k = m^k revisited using continued fractions,
// Math. Comp. 80 (2011), 1221-1237.

class CF
{
private:
	Heap & _heap;
	guint _N, _cond_b, _a_j;
	uint64_t _j;
	gfloat _q_j, _q_jm1;
	uint64_t _q_j_mod, _q_jm1_mod;
	std::string _Nstr;
	Factor _factor;
	// { p : 3 is a primitive root modulo p } such that 6 * (5*7*17*19*29*31*43)^2 < 2^64
	const std::vector<uint32_t> _P_pr3 = { 5, 7, 17, 19, 29, 31, 43 };
	static const uint64_t _mod_q = 6 * uint64_t(5*7*17*19*29*31*43) * (5*7*17*19*29*31*43);
	// uint64_t mod_q_inv; int mod_q_e; gint::mod_init(_mod_q, mod_q_inv, mod_q_e);
	static const uint64_t _mod_q_inv = 2319961602461247110ull;
	static const int _mod_q_e = 60;
	static const uint64_t _mod_q_f = (-_mod_q) % _mod_q;	// 2^64 mod mod_q

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

	void _gcf_extend(Mat22u & M, guint & d, const uint64_t n, const uint64_t size)
	{
		if (size == 64)
		{
			M.set_gcf(n);
			for (uint64_t i = 1; i < size; ++i) M.mul_right_gcf(n + i);

			d = 1;
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
			_gcf_extend(M, d, n, size / 2);
			Mat22u M_r; guint d_r; _gcf_extend(M_r, d_r, n + size / 2, size / 2);
			M.mul_right(M_r); d *= d_r;
		}
	}

	double gcf_extend(Mat22u & M, guint & divisor, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_gcf_extend(M, divisor, n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	static double gcf_mul(Mat22 & M, const Mat22u & Mgcf)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.mul_right(Mgcf);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	static double gcf_div(Mat22 & M, guint & divisor)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.div(divisor);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	size_t _cf_reduce(Mat22 & M, Mat22 & Mcf, bool & found)
	{
		size_t count = 0;
		Mcf.set_identity();

		// a_j is the jth coefficient of the regular continued fraction
		guint a_jp1, a_jp2;
		while (M.get_cf_coefficient(a_jp1, a_jp2))
		{
			Mcf.cf_mul(a_jp1); Mcf.cf_mul(a_jp2);

			_j += 2;	// j must always be odd, now a_{j+2} is a_j
			++count;

			if (a_jp2.cmp(_cond_b) >= 0)	// We have: (a): j-1 is even, (b) a_j >= 180N - 2
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
		const size_t n = M.get_min_size();
		if (n < 32) return _cf_reduce(M, Mcf, found);

		Mat22 M_lo; M.split(M_lo, n / 2);

		const size_t count1 = _cf_reduce_half(level + 1, M, Mcf, found);
		if (count1 == 0) { M.combine(M_lo, Mcf, n / 2); return 0; }

		// M_hi is Mcf * M.hi then Mcf * M = (M_hi << (n/2 * GMP_LIMB_BITS)) + Mcf * M_lo
		M_lo.mul_left(Mcf); M.combine(M_lo, Mcf, n / 2);

		if (found) return count1;

		const size_t n2 = M.get_min_size();
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
		const uint64_t i11 = Mcf.get11().mod(_mod_q, _mod_q_inv, _mod_q_e, _mod_q_f), i12 = Mcf.get12().mod(_mod_q, _mod_q_inv, _mod_q_e, _mod_q_f);
		const uint64_t i21 = Mcf.get21().mod(_mod_q, _mod_q_inv, _mod_q_e, _mod_q_f), i22 = Mcf.get22().mod(_mod_q, _mod_q_inv, _mod_q_e, _mod_q_f);
		Mod64 mod64(_mod_q);
		const uint64_t ti = mod64.sub(mod64.mul(i11, _q_jm1_mod), mod64.mul(i12, _q_j_mod));
		_q_j_mod = mod64.sub(mod64.mul(i22, _q_j_mod), mod64.mul(i21, _q_jm1_mod)); _q_jm1_mod = ti;

		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool condition_c() const
	{
		const uint64_t q_j_mod = _q_jm1_mod;	// a step backward
		const uint64_t g6 = q_j_mod % 6;
		return ((g6 == 1) || (g6 == 5));	// (c) gcd(q_{j , 6) = 1
	}

	bool condition_d() const
	{
		// a step backward
		const uint64_t j = _j - 1;
		const guint & a_jp1 = _a_j;
		const gfloat & q_j = _q_jm1;
		const uint64_t q_j_mod = _q_jm1_mod;

		const uint64_t g6 = q_j_mod % 6;

		std::ostringstream ss, ssr;
		ss << "N = " << _Nstr << ", j = " << j << ", a_{j+1} = " << a_jp1.to_string() << ", q_j = " << q_j.to_string()
			<< ", q_j mod 6 = " << ((g6 == 1) ? "+" : "-") << "1";
		ssr << _Nstr << "\t" << j << "\t" << a_jp1.to_string() << "\t" << q_j.to_string() << "\t" << ((g6 == 1) ? "+" : "-") << "1";

		bool cond_d = true;
		for (const uint32_t p : _P_pr3)
		{
			const uint64_t r = q_j_mod % (p * p);
			if ((r != 0) && (r % p == 0)) { cond_d = false; ss << ", (d) p = " << p; ssr << "\t" << p; break; }
		}
		ss << "." << std::endl; ssr << std::endl;
		record(ss.str(), ssr.str());

		return cond_d;
	}

	std::string checkpoint_filename() const { return std::string("cf_") + _N.to_string() + ".cpt"; }

	int _read_checkpoint(const std::string & filename, uint64_t & n, uint64_t & nstep, Mat22 & M)
	{
		Checkpoint checkpoint(filename, "rb");
		if (!checkpoint.exists()) return -1;

		uint32_t version = 0;
		if (!checkpoint.read(reinterpret_cast<char *>(&version), sizeof(version))) return -2;
		if (version != 1) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&_j), sizeof(_j))) return -2;
		double mantissa; size_t exponent;
		if (!checkpoint.read(reinterpret_cast<char *>(&mantissa), sizeof(mantissa))) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&exponent), sizeof(exponent))) return -2;
		_q_j = gfloat(mantissa, exponent);
		if (!checkpoint.read(reinterpret_cast<char *>(&mantissa), sizeof(mantissa))) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&exponent), sizeof(exponent))) return -2;
		_q_jm1 = gfloat(mantissa, exponent);
		if (!checkpoint.read(reinterpret_cast<char *>(&_q_j_mod), sizeof(_q_j_mod))) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&_q_jm1_mod), sizeof(_q_jm1_mod))) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&n), sizeof(n))) return -2;
		if (!checkpoint.read(reinterpret_cast<char *>(&nstep), sizeof(nstep))) return -2;
		if (!M.read(checkpoint)) return -2;
		if (!checkpoint.check_crc32()) return -2;
		return 0;
	}

	bool read_checkpoint(uint64_t & n, uint64_t & nstep, Mat22 & M)
	{
		std::string cpt_file = checkpoint_filename();
		int error = _read_checkpoint(cpt_file, n, nstep, M);
		if (error < -1) std::cerr << "Error: " << cpt_file << " is an invalid checkpoint." << std::endl;
		if (error < 0)
		{
 			cpt_file += ".old";
			error = _read_checkpoint(cpt_file, n, nstep, M);
			if (error < -1) std::cerr << "Error: " << cpt_file << " is an invalid checkpoint." << std::endl;
		}
		return (error == 0);
	}

	void save_checkpoint(const uint64_t n, const uint64_t nstep, const Mat22 & M)
	{
		const std::string cpt_file = checkpoint_filename(), old_cpt_file = cpt_file + ".old", new_cpt_file = cpt_file + ".new";

		{
			Checkpoint checkpoint(new_cpt_file, "wb");
			uint32_t version = 1;
			if (!checkpoint.write(reinterpret_cast<const char *>(&version), sizeof(version))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&_j), sizeof(_j))) return;
			const double q_j_mantissa = _q_j.get_mantissa(); const size_t q_j_exponent = _q_j.get_exponent();
			if (!checkpoint.write(reinterpret_cast<const char *>(&q_j_mantissa), sizeof(q_j_mantissa))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&q_j_exponent), sizeof(q_j_exponent))) return;
			const double _q_jm1_mantissa = _q_jm1.get_mantissa(); const size_t _q_jm1_exponent = _q_jm1.get_exponent();
			if (!checkpoint.write(reinterpret_cast<const char *>(&_q_jm1_mantissa), sizeof(_q_jm1_mantissa))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&_q_jm1_exponent), sizeof(_q_jm1_exponent))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&_q_j_mod), sizeof(_q_j_mod))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&_q_jm1_mod), sizeof(_q_jm1_mod))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&n), sizeof(n))) return;
			if (!checkpoint.write(reinterpret_cast<const char *>(&nstep), sizeof(nstep))) return;
			if (!M.write(checkpoint)) return;
			checkpoint.write_crc32();
		}

		std::remove(old_cpt_file.c_str());

		struct stat s;
		if ((stat(cpt_file.c_str(), &s) == 0) && (std::rename(cpt_file.c_str(), old_cpt_file.c_str()) != 0))	// file exists and cannot rename it
		{
			std::cerr << "Error: cannot save checkpoint." << std::endl;
			return;
		}

		if (std::rename(new_cpt_file.c_str(), cpt_file.c_str()) != 0)
		{
			std::cerr << "Error: cannot save checkpoint." << std::endl;
			return;
		}
	}

	void record(const std::string & str, const std::string & str_res = "") const
	{
		std::cout << str;
		std::ofstream logfile("cflog.txt", std::ios::app);
		if (logfile.is_open())
		{
			logfile << str;
			logfile.close();
		}
		if (!str_res.empty())
		{
			std::ofstream resfile("cfres.txt", std::ios::app);
			if (resfile.is_open())
			{
				resfile << str_res;
				resfile.close();
			}
		}
	}

public:
	void solve(const guint & N, const std::string & Nstr, const bool verbose)
	{
		_heap.reset();

		_N = N;
		_cond_b = N; _cond_b *= 180; _cond_b -= 2;
		_Nstr = Nstr;

		_factor.init();

		// n is the index of the convergent of the generalized continued fraction.
		uint64_t n, nstep;

		// Matrix M is the remainder of the nth convergent the generalized continued fraction of log(2) / 2N
		// after j coefficients of the regular continued fraction retrieval using the Euclidean algorithm.
		// j and n must be odd such that a_11/a21 < a12/a22
		Mat22 M;

		const bool resume = read_checkpoint(n, nstep, M);
		if (resume)
		{
			std::cout << "Resuming from a checkpoint." << std::endl;
		}
		else
		{
			// j is the index of the convergent of the regular continued fraction.
			_j = uint64_t(-1);

			// Denominator of the regular continued fraction: q_j and q_{j-1}.
			_q_j = gfloat(0, 0); _q_jm1 = gfloat(1, 0);
			_q_j_mod = 0; _q_jm1_mod = 1;

			n = 1; nstep = 256;

			// p_0 = 0, p_1 = 1, q_0 = q_1 = 2N.
			M.init_gcf(N);
		}

		double time_gcf_extend = 0, time_gcf_mul = 0, time_gcf_div = 0, time_cf_reduce = 0;
		size_t M_min_size = 0, Mgcf_size = 0;	// divisor_size = 0, M_max_size = 0;
		uint64_t j_prev = _j;
		size_t str_size = 0;

		const auto start_time = std::chrono::high_resolution_clock::now();
		auto display_time = start_time, save_time = start_time;

		bool found = false;
		while (!found)
		{
			{
				// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
				// Compute nstep terms of the matrix and its divisor starting at n.
				Mat22u Mgcf; guint divisor;
				time_gcf_extend += gcf_extend(Mgcf, divisor, n, nstep);
				Mgcf_size = Mgcf.get_byte_count();
				// divisor_size = divisor.get_byte_count();

				time_gcf_mul += gcf_mul(M, Mgcf);
				time_gcf_div += gcf_div(M, divisor);
			}
			n += nstep;

			// M_max_size = M.get_byte_count();

			time_cf_reduce += cf_reduce(M, found);

			M_min_size = M.get_byte_count();

			if (Mgcf_size < M_min_size / 4) nstep *= 2;

			if (found) found = condition_c();

			const auto now = std::chrono::high_resolution_clock::now();

			if (std::chrono::duration<double>(now - save_time).count() > 10) { save_time = now; save_checkpoint(n, nstep, M); }

			if ((std::chrono::duration<double>(now - display_time).count() > 3600) || found)
			{
				display_time = now;
				const double elapsed_time = std::chrono::duration<double>(now - start_time).count();

				std::ostringstream ss;
				ss	<< "N = " << Nstr << ", n = " << n << " [" << nstep << "], j = " << _j << " (+" << _j - j_prev << "), " << "q_j = " << _q_j.to_string();
				if (!verbose) ss << ", memory usage: " << _heap.get_memory_usage();
				ss << ", elapsed time: " << format_time(elapsed_time) << ".";

				if (verbose)
				{
					const double time_total = time_gcf_extend + time_gcf_mul + time_gcf_div + time_cf_reduce;

					ss << std::endl;
					// << "M_max: " << M_max_size << ", M_min: " << M_min_size << ", Mgcf: " << Mgcf_size << ", divisor: " << divisor_size << ", "
					ss	<< "    Memory usage: " << _heap.get_memory_info() << std::endl;
					ss	<< "    CPU usage: " << std::setprecision(3)
						<< "gcf_extend: " << time_gcf_extend * 100 / time_total << "%, gcf_mul: " << time_gcf_mul * 100 / time_total << "%, "
						<< "gcf_div: " << time_gcf_div * 100 / time_total << "%, cf_reduce: " << time_cf_reduce * 100 / time_total << "%." << std::endl;
				}
				else
				{
					const size_t s_size = ss.str().size();
					for (size_t i = s_size; i < str_size; ++i) ss << " ";
					str_size = s_size;
					if (found) ss << std::endl; else ss << "\r";
				}

				record(ss.str());
				j_prev = _j;
			}

			if (found) found = condition_d();
		}

		_N = 0; _cond_b = 0; _a_j = 0; M.set_zero();
		std::cout << "Memory size: " << _heap.get_memory_size() << "." << std::endl;
	}
};
