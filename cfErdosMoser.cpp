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

class Mat22
{
private:
	gint _a11, _a12;
	gint _a21, _a22;

public:
	Mat22() {}
	virtual ~Mat22() {}
	Mat22(const Mat22 & rhs) : _a11(rhs._a11), _a12(rhs._a12), _a21(rhs._a21), _a22(rhs._a22) {}

	Mat22 & operator = (const Mat22 & rhs)
	{
		if (&rhs == this) return *this;
		_a11 = rhs._a11; _a12 = rhs._a12; _a21 = rhs._a21; _a22 = rhs._a22;
		return *this;
	}

	const gint & get11() const { return _a11; }
	const gint & get12() const { return _a12; }
	const gint & get21() const { return _a21; }
	const gint & get22() const { return _a22; }

	void set_identity() { _a11 = 1u; _a12 = 0u; _a21 = 0u; _a22 = 1u; }

	void set_gcf(const uint64_t n)
	{
		_a11 = n; _a12 = 2 * n + 1; _a12 *= _a11;
		_a21 = 2u; _a22 = 5 * n + 2;
	}

private:
	void _mul_right(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a12);
		_a11 *= rhs._a11; _a11.addmul(_a12, rhs._a21);
		t.addmul(_a12, rhs._a22); _a12.swap(t);

		t.mul(_a21, rhs._a12);
		_a21 *= rhs._a11; _a21.addmul(_a22, rhs._a21);
		t.addmul(_a22, rhs._a22); _a22.swap(t);
	}

	void _mul_left(const Mat22 & rhs)
	{
		gint t;

		t.mul(_a11, rhs._a21);
		_a11 *= rhs._a11; _a11.addmul(_a21, rhs._a12);
		t.addmul(_a21, rhs._a22); _a21.swap(t);

		t.mul(_a12, rhs._a21);
		_a12 *= rhs._a11; _a12.addmul(_a22, rhs._a12);
		t.addmul(_a22, rhs._a22); _a22.swap(t);
	}

	void _div(const gint & d) { _a11.divexact(d); _a12.divexact(d); _a21.divexact(d); _a22.divexact(d); }

public:
	size_t get_byte_count() const { return _a11.get_byte_count() + _a12.get_byte_count() + _a21.get_byte_count() + _a22.get_byte_count(); }

	double mul_right(const Mat22 & rhs)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_mul_right(rhs);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void mul_left(const Mat22 & rhs)
	{
		_mul_left(rhs);
	}

	double div(const gint & d)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_div(d);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void split(const Mat22 & rhs, const size_t n)
	{
		_a11.rshift(rhs._a11, n);
		_a12.rshift(rhs._a12, n); _a12 += 1;
		_a21.rshift(rhs._a21, n); _a21 += 1;
		_a22.rshift(rhs._a22, n);
	}

	void init_gcf(const gint & N)
	{
		_a11 = 0u; _a12 = 1u;
		_a21 = N; _a21 += _a21; _a22 = _a21;
	}

	void cf_mul(const gint & coefficient)
	{
		_a11.submul(coefficient, _a21); _a11.swap(_a21);
		_a12.submul(coefficient, _a22); _a12.swap(_a22);
	}

	bool get_cf_coefficient(gint & coefficient)
	{
		coefficient.div(_a11, _a21);
		gint t; t.div(_a12, _a22);
		const bool success = (coefficient == t);
		if (success) cf_mul(coefficient);
		return success;
	}
};

class CF
{
private:
	gint _N, _cond_b, _q_j, _q_jm1, _a_j, _t;
	std::vector<uint32_t> _PN;
	uint64_t _j;
	Factor factor;

public:
	CF() {}
	virtual ~CF() {}

private:
	void _gcf_matrix(Mat22 & M, const uint64_t n, const uint64_t size)
	{
		if (size == 8)
		{
			M.set_gcf(n);
			Mat22 Mi;
			for (uint64_t i = 1; i < size; ++i)
			{
				Mi.set_gcf(n + i);
				M.mul_right(Mi);
			}
		}
		else
		{
			_gcf_matrix(M, n, size / 2);
			Mat22 Mr; _gcf_matrix(Mr, n + size / 2, size / 2);
			M.mul_right(Mr);
		}
	}

	double gcf_matrix(Mat22 & M, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_gcf_matrix(M, n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void _cfg_divisor(gint & divisor, const uint64_t n, const uint64_t size)
	{
		if (size == 8)
		{
			divisor = 1u;
			for (uint64_t i = 0; i < size; ++i)
			{
				const uint64_t n_i = n + i;
				const uint64_t p = factor.smallest(n_i);
				// if n = p^k is a prime power then pp = p else pp = 1 (exponential of Mangoldt function).
				uint64_t m = 1; if (p != n_i) { m = n_i; while (m % p == 0) m /= p; }
				uint64_t pp = (m == 1) ? p : 1;
				_t = n_i / pp;
				divisor *= _t;
			}
		}
		else
		{
			_cfg_divisor(divisor, n, size / 2);
			gint divisor_r; _cfg_divisor(divisor_r, n + size / 2, size / 2);
			divisor *= divisor_r;
		}
	}

	double cfg_divisor(gint & divisor, const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_cfg_divisor(divisor, n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	double gcf_mul(Mat22 & M, const Mat22 & Mg)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.mul_right(Mg);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	double gcf_div(Mat22 & M, const gint & divisor)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		M.div(divisor);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool _direct(const Mat22 & uv, Mat22 & M, bool & found)
	{
		bool succeed = false;

		M.set_identity();
		Mat22 uvp(uv);

		// a_j is the jth coefficient of the regular continued fraction
		while (uvp.get_cf_coefficient(_a_j))
		{
			// Update the regular continued fraction
			_q_jm1.addmul(_a_j, _q_j); _q_j.swap(_q_jm1);

			M.cf_mul(_a_j);

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
							const uint32_t nu_q_jm1 = _t.nu(_q_jm1, p);
							if (nu_q_jm1 > 0)	// p | q_{j - 1}
							{
								// nu(3^{p−1} − 1) = 1 because p != 11, 1006003 (Mirimanoff primes)
								if (nu_q_jm1 != 1 + _t.nu(_N, p) + 1) { cond_d = false; std::cout << ", " << p; }	// (d)
							}
						}

						std::cout << std::endl;

						if (cond_d) { found = true; break; }
					}
				}
			}

			++_j;
			succeed = true;
		}

		return succeed;
	}

	bool _half(const Mat22 & uv, Mat22 & M, bool & found)
	{
		const size_t n = uv.get12().get_word_count();

		if (n <= 32)
		{
			return _direct(uv, M, found);
		}

		Mat22 uv_h; uv_h.split(uv, n / 2);

		if (!_half(uv_h, M, found)) return false;
		if (found) return true;

		uv_h = uv;
		uv_h.mul_left(M);
		uv_h.split(uv_h, n / 2);

		Mat22 Mp;
		if (_half(uv_h, Mp, found))
		{
			if (!found) M.mul_left(Mp);
		}

		return true;
	}

	double half(Mat22 & uv, bool & found)
	{
		const auto start = std::chrono::high_resolution_clock::now();

		const size_t n = uv.get12().get_word_count();	// bytes

		if (n <= 32)
		{
			Mat22 M; _direct(uv, M, found);
			uv.mul_left(M);
		}
		else
		{
			Mat22 nuv;
			while (true)
			{
				nuv.split(uv, n / 2);

				Mat22 M;
				if (!_half(nuv, M, found)) break;
				if (found) break;

				uv.mul_left(M);
			}
			// nuv.split(uv, n / 2);

			// if (_half(nuv, M, found))
			// {
			// 	if (!found) uv.mul_left(M);
			// }
		}

		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

public:
	void solve(const gint & N, const std::vector<uint32_t> & PN)
	{
		_N = N; _PN = PN;
		_cond_b = N; _cond_b *= 180; _cond_b -= 2u;

		factor.init();

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

		double time_gcf_matrix = 0, time_gcf_divisor = 0, time_gcf_mul = 0, time_gcf_div = 0, time_mul_left = 0;
		size_t M_size = 0, Mg_size = 0;

		bool found = false;
		while (!found)
		{
			{
				// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
				// Compute nstep terms starting at n.
				Mat22 Mg; time_gcf_matrix += gcf_matrix(Mg, n, nstep);
				// divisor of M
				gint divisor; time_gcf_divisor += cfg_divisor(divisor, n, nstep);

				time_gcf_mul += gcf_mul(M, Mg);
				time_gcf_div += gcf_div(M, divisor);

				Mg_size = std::max(Mg_size, Mg.get_byte_count());
			}
			n += nstep;

			M_size = std::max(M_size, M.get_byte_count()); 

			time_mul_left += half(M, found);

			if (Mg_size < M.get_byte_count() / 2) nstep *= 2;

			const double time_elapsed = time_gcf_matrix + time_gcf_divisor + time_gcf_mul + time_gcf_div + time_mul_left;
			std::cout << std::setprecision(3)
				<< "j: " << _j << ", elapsed time: " << time_elapsed << ", "
				<< "M_max: " << M_size << ", M_min: " << M.get_byte_count() << ", Mg: " << Mg_size << ", q_j: " << _q_j.get_byte_count() << ", "
				<< "gcf_matrix: " << time_gcf_matrix * 100 / time_elapsed << "%, "
				<< "gcf_divisor: " << time_gcf_divisor * 100 / time_elapsed << "%, "
				<< "gcf_mul: " << time_gcf_mul * 100 / time_elapsed << "%, "
				<< "gcf_div: " << time_gcf_div * 100 / time_elapsed << "%, "
				<< "mul_left: " << time_mul_left * 100 / time_elapsed << "%."
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
