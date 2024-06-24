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

#include <gmp.h>

inline void mpz_set_ui_64(mpz_t rop, const uint64_t n)
{
#ifdef _LONG_LONG_LIMB
	mpz_set_ui(rop, (unsigned long int)n); if (n > 0) rop->_mp_d[0] = n;
#else
	mpz_set_ui(rop, (unsigned long int)(n >> 32));
	mpz_mul_2exp(rop, rop, 32);
	mpz_add_ui(rop, rop, (unsigned long int)n);
#endif
}

#include "factor.h"

class Mat22
{
private:
	mpz_t _a11, _a12;
	mpz_t _a21, _a22;

public:
	Mat22() { mpz_inits(_a11, _a12, _a21, _a22, nullptr); }
	virtual ~Mat22() { mpz_clears(_a11, _a12, _a21, _a22, nullptr); }
	Mat22(const Mat22 & rhs) { mpz_init_set(_a11, rhs._a11); mpz_init_set(_a12, rhs._a12); mpz_init_set(_a21, rhs._a21); mpz_init_set(_a22, rhs._a22); }

	Mat22 & operator = (const Mat22 & rhs)
	{
		if (&rhs == this) return *this;
		mpz_set(_a11, rhs._a11); mpz_set(_a12, rhs._a12); mpz_set(_a21, rhs._a21); mpz_set(_a22, rhs._a22);
		return *this;
	}

	const mpz_t & get11() const { return _a11; }
	const mpz_t & get12() const { return _a12; }
	const mpz_t & get21() const { return _a21; }
	const mpz_t & get22() const { return _a22; }

	void set_identity() { mpz_set_ui(_a11, 1); mpz_set_ui(_a12, 0); mpz_set_ui(_a21, 0); mpz_set_ui(_a22, 1); }

private:
	void _mul_right(const Mat22 & rhs)
	{
		mpz_t t; mpz_init(t);

		mpz_mul(t, _a11, rhs._a12);
		mpz_mul(_a11, _a11, rhs._a11); mpz_addmul(_a11, _a12, rhs._a21);
		mpz_addmul(t, _a12, rhs._a22); mpz_swap(_a12, t);

		mpz_mul(t, _a21, rhs._a12);
		mpz_mul(_a21, _a21, rhs._a11); mpz_addmul(_a21, _a22, rhs._a21);
		mpz_addmul(t, _a22, rhs._a22); mpz_swap(_a22, t);

		mpz_clear(t);
	}

	void _mul_left(const Mat22 & rhs)
	{
		mpz_t t; mpz_init(t);

		mpz_mul(t, _a11, rhs._a21);
		mpz_mul(_a11, _a11, rhs._a11); mpz_addmul(_a11, _a21, rhs._a12);
		mpz_addmul(t, _a21, rhs._a22); mpz_swap(_a21, t);

		mpz_mul(t, _a12, rhs._a21);
		mpz_mul(_a12, _a12, rhs._a11); mpz_addmul(_a12, _a22, rhs._a12);
		mpz_addmul(t, _a22, rhs._a22); mpz_swap(_a22, t);

		mpz_clear(t);
	}

	void _div(const mpz_t & d)
	{
		static const std::string error = "divisor";
		if (mpz_divisible_p(_a11, d)) mpz_divexact(_a11, _a11, d); else throw std::runtime_error(error);
		if (mpz_divisible_p(_a12, d)) mpz_divexact(_a12, _a12, d); else throw std::runtime_error(error);
		if (mpz_divisible_p(_a21, d)) mpz_divexact(_a21, _a21, d); else throw std::runtime_error(error);
		if (mpz_divisible_p(_a22, d)) mpz_divexact(_a22, _a22, d); else throw std::runtime_error(error);
	}

	void _eval_gcf(const uint64_t n, const uint64_t size)
	{
		if (size == 1)
		{
			mpz_set_ui_64(_a11, n); mpz_set_ui_64(_a12, 2 * n + 1); mpz_mul(_a12, _a12, _a11);
			mpz_set_ui(_a21, 2); mpz_set_ui_64(_a22, 5 * n + 2);
		}
		else if (size == 8)
		{
			_eval_gcf(n, 1);
			Mat22 M;
			for (uint64_t i = 1; i < size; ++i)
			{
				M._eval_gcf(n + i, 1);
				mul_right(M);
			}
		}
		else
		{
			_eval_gcf(n, size / 2);
			Mat22 M; M._eval_gcf(n + size / 2, size / 2);
			mul_right(M);
		}
	}

public:
	size_t size() const { return (mpz_size(_a11) + mpz_size(_a12) + mpz_size(_a21) + mpz_size(_a22)) * sizeof(mp_limb_t); }

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

	double div(const mpz_t & d)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_div(d);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void split(const Mat22 & rhs, const size_t n)
	{
		mpz_div_2exp(_a11, rhs._a11, n * sizeof(mp_limb_t));
		mpz_div_2exp(_a12, rhs._a12, n * sizeof(mp_limb_t)); mpz_add_ui(_a12, _a12, 1);
		mpz_div_2exp(_a21, rhs._a21, n * sizeof(mp_limb_t)); mpz_add_ui(_a21, _a21, 1);
		mpz_div_2exp(_a22, rhs._a22, n * sizeof(mp_limb_t));
	}

	void init_gcf(const mpz_t & N)
	{
		mpz_set_ui(_a11, 0); mpz_set_ui(_a12, 1);
		mpz_set(_a21, N); mpz_mul_2exp(_a21, _a21, 1); mpz_set(_a22, _a21);
	}

	double eval_gcf(const uint64_t n, const uint64_t size)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_eval_gcf(n, size);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	void cf_mul(const mpz_t & coefficient)
	{
		mpz_submul(_a11, coefficient, _a21); mpz_swap(_a11, _a21);
		mpz_submul(_a12, coefficient, _a22); mpz_swap(_a12, _a22);
	}

	bool get_cf_coefficient(mpz_t & coefficient)
	{
		mpz_t t; mpz_init(t);
		mpz_fdiv_q(coefficient, _a11, _a21);
		mpz_fdiv_q(t, _a12, _a22);
		const bool success = (mpz_cmp(coefficient, t) == 0);
		mpz_clear(t);
		if (success) cf_mul(coefficient);
		return success;
	}
};

class CF
{
private:
	mpz_t _N, _cond_b, _q_j, _q_jm1, _mdivisor, _a_j, _t;
	std::vector<uint32_t> _PN;
	uint64_t _j;
	Factor factor;

public:
	CF() { mpz_inits(_N, _cond_b, _q_j, _q_jm1, _mdivisor, _a_j, _t, nullptr); }
	virtual ~CF() { mpz_clears(_N, _cond_b, _q_j, _q_jm1, _mdivisor, _a_j, _t, nullptr); }

private:
	uint32_t nu(const mpz_t & n, const uint32_t p)
	{
		mpz_set(_t, n);
		uint32_t a = 0; while (mpz_divisible_ui_p(_t, p)) { mpz_divexact_ui(_t, _t, p); ++a; }
		return a;
	}

	void _eval_mdivisor(const uint64_t n, const uint64_t size, mpz_t & mdivisor)
	{
		if (size == 8)
		{
			mpz_set_ui(mdivisor, 1);
			for (uint64_t i = 0; i < size; ++i)
			{
				const uint64_t n_i = n + i;
				const uint64_t p = factor.smallest(n_i);
				// if n = p^k is a prime power then pp = p else pp = 1 (exponential of Mangoldt function).
				uint64_t m = 1; if (p != n_i) { m = n_i; while (m % p == 0) m /= p; }
				uint64_t pp = (m == 1) ? p : 1;
				mpz_set_ui_64(_t, n_i / pp);
				mpz_mul(mdivisor, mdivisor, _t);
			}
		}
		else
		{
			_eval_mdivisor(n, size / 2, mdivisor);
			mpz_t rdivisor; mpz_init(rdivisor); _eval_mdivisor(n + size / 2, size / 2, rdivisor);
			mpz_mul(mdivisor, mdivisor, rdivisor);
			mpz_clear(rdivisor);
		}
	}

	double eval_mdivisor(const uint64_t n, const uint64_t size, mpz_t & mdivisor)
	{
		const auto start = std::chrono::high_resolution_clock::now();
		_eval_mdivisor(n, size, mdivisor);
		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

	bool direct(const Mat22 & uv, Mat22 & M, bool & found)
	{
		bool succeed = false;

		M.set_identity();
		Mat22 uvp(uv);

		// a_j is the jth coefficient of the regular continued fraction
		while (uvp.get_cf_coefficient(_a_j))
		{
			// Update the regular continued fraction
			mpz_addmul(_q_jm1, _a_j, _q_j); mpz_swap(_q_j, _q_jm1);

			M.cf_mul(_a_j);

			if (_j % 2 == 1)	// (a): j - 1 is even
			{
				if (mpz_cmp(_a_j, _cond_b) >= 0)	// (b) a_j >= 180N - 2
				{
					const unsigned long int g6 = mpz_fdiv_ui(_q_jm1, 6);
					if ((g6 == 1) || (g6 == 5))		// (c) gcd(q_{j - 1} , 6) = 1
					{
						// time_mul_left += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();

						long e2; const double x = mpz_get_d_2exp(&e2, _q_jm1);
						const double lx = std::log10(x) + e2 * std::log10(2.0);
						const uint64_t e10 = uint64_t(lx); const double m10 = std::pow(10.0, lx - double(e10));

						mpz_out_str(stdout, 10, _N);
						std::cout << ", " << _j - 1 << ", ";
						mpz_out_str(stdout, 10, _a_j);
						std::cout << ", " << m10 << "*10^" << e10;
						std::cout << ", " << ((g6 == 1) ? "+" : "-") << "1";

						bool cond_d = true;
						for (const uint32_t p : _PN)
						{
							const uint32_t nu_q_jm1 = nu(_q_jm1, p);
							if (nu_q_jm1 > 0)	// p | q_{j - 1}
							{
								// nu(3^{p−1} − 1) = 1 because p != 11, 1006003 (Mirimanoff primes)
								if (nu_q_jm1 != 1 + nu(_N, p) + 1) { cond_d = false; std::cout << ", " << p; }	// (d)
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
		const size_t n = mpz_size(uv.get12());

		if (n <= 32)
		{
			return direct(uv, M, found);
		}

		Mat22 nuv;
		nuv.split(uv, n / 2);

		const bool succeed = _half(nuv, M, found);
		if (!succeed) return false;
		if (found) return true;

		nuv = uv;
		nuv.mul_left(M);
		nuv.split(nuv, n / 2);

		Mat22 Mp;
		if (_half(nuv, Mp, found))
		{
			if (found) return true;
			M.mul_left(Mp);
		}

		return true;
	}

	double half(Mat22 & uv, bool & found)
	{
		const auto start = std::chrono::high_resolution_clock::now();

		const size_t n = mpz_size(uv.get12());

		if (n <= 32)
		{
			Mat22 M; direct(uv, M, found);
			uv.mul_left(M);
		}
		else
		{
			Mat22 nuv;
			nuv.split(uv, n / 2);

			Mat22 M;
			if (_half(nuv, M, found))
			{
				if (!found)
				{
					uv.mul_left(M);

					nuv.split(uv, n / 2);

					if (_half(nuv, M, found))
					{
						if (!found) uv.mul_left(M);
					}
				}
			}
		}

		return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
	}

public:
	void solve(const mpz_t & N, const std::vector<uint32_t> & PN)
	{
		mpz_set(_N, N); _PN = PN;
		mpz_mul_ui(_cond_b, N, 180); mpz_sub_ui(_cond_b, _cond_b, 2);

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
		mpz_set_ui(_q_j, 0); mpz_set_ui(_q_jm1, 1);

		double time_eval_gcf = 0, time_eval_mdivisor = 0, time_mul_right = 0, time_div = 0, time_mul_left = 0;
		size_t M_size = 0, Mg_size = 0;

		bool found = false;
		while (!found)
		{
			// Matrix form of the generalized continued fraction: p_{n-2} / q_{n-2} and p_{n-1} / q_{n-1}.
			// Compute nstep terms starting at n.
			Mat22 Mg; time_eval_gcf += Mg.eval_gcf(n, nstep);
			time_eval_mdivisor += eval_mdivisor(n, nstep, _mdivisor);
			n += nstep;

			time_mul_right += M.mul_right(Mg);
			time_div += M.div(_mdivisor);

			M_size = std::max(M_size, M.size()); Mg_size = std::max(Mg_size, Mg.size());

			time_mul_left += half(M, found);

			if (Mg_size < M.size() / 2) nstep *= 2;

			const double time_elapsed = time_eval_gcf + time_eval_mdivisor + time_mul_right + time_div + time_mul_left;
			std::cout << std::setprecision(3)
				<< "j: " << _j << ", elapsed time: " << time_elapsed << ", "
				<< "M_max: " << M_size << ", M_min: " << M.size() << ", Mg: " << Mg_size << ", q_j: " << mpz_size(_q_j) * sizeof(mp_limb_t) << ", "
				<< "eval_gcf: " << time_eval_gcf * 100 / time_elapsed << "%, "
				<< "eval_mdivisor: " << time_eval_mdivisor * 100 / time_elapsed << "%, "
				<< "mul_right: " << time_mul_right * 100 / time_elapsed << "%, "
				<< "div: " << time_div * 100 / time_elapsed << "%, "
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
		mpz_t N; mpz_init(N);

		// for (size_t d = 0; d < fN.size(); ++d)
		size_t d = 10;
		{
			mpz_set_ui(N, 1);
			for (size_t i = 0; i < d; ++i) mpz_mul_ui(N, N, fN[i]);

			std::vector<uint32_t> PN;	// PN is a subset of PNmax
			for (const uint32_t p : PNmax) if (mpz_divisible_ui_p(N, p - 1)) PN.push_back(p);

			cf.solve(N, PN);
		}

		mpz_clear(N);
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << "." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
