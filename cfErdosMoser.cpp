/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>
#include <set>

#ifdef _WIN32
#include <Windows.h>
#else
#include <signal.h>
#endif

#include "CF.h"

#include <gmp.h>

class Application
{
private:
	struct deleter { void operator()(const Application * const p) { delete p; } };

	static void quit(int) { CF::get_instance().quit(); }

#ifdef _WIN32
	static BOOL WINAPI HandlerRoutine(DWORD) { quit(1); return TRUE; }
#endif

public:
	Application()
	{
#ifdef _WIN32
		SetConsoleCtrlHandler(HandlerRoutine, TRUE);
#else
		signal(SIGTERM, quit);
		signal(SIGINT, quit);
#endif
	}

	virtual ~Application() {}

	static Application & get_instance()
	{
		static std::unique_ptr<Application, deleter> p_instance(new Application());
		return *p_instance;
	}

private:
	static std::string header()
	{
		const char * const sysver =
#if defined(_WIN64)
			"win64";
#elif defined(_WIN32)
			"win32";
#elif defined(__linux__)
#ifdef __x86_64
			"linux64";
#else
			"linux32";
#endif
#elif defined(__APPLE__)
			"macOS";
#else
			"unknown";
#endif

		std::ostringstream ssc;
#if defined(__GNUC__)
		ssc << " gcc-" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
#elif defined(__clang__)
		ssc << " clang-" << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#endif

		std::ostringstream ss;
		ss << "cfErdosMoser 24.07.0 " << sysver << ssc.str() << std::endl;
		ss << "Copyright (c) 2024, Yves Gallot" << std::endl;
		ss << "cfErdosMoser is free source code, under the MIT license." << std::endl << std::endl;
		return ss.str();
	}

	static std::string usage()
	{
		std::ostringstream ss;
		ss << "Usage: cfErdosMoser <N> [-v]" << std::endl;
		ss << "  N: compute the regular continued fraction of log(2)/(2N)," << std::endl;
		ss << "     if N = 0 then check N = 1, 2, 2^2, ..., 2^8, 2^8*3, ..." << std::endl;
		ss << " -t <n>: number of threads (default: 4 threads, 0: all logical cores)," << std::endl;
		ss << " -v: verbose mode." << std::endl;
		return ss.str();
	}

public:
	void run(int argc, char * argv[])
	{
		// for (size_t z_size = 256; z_size <= (size_t(1) << 32); z_size *= 2)
		// {
		// 	int e = ((z_size & (~z_size + 1)) == z_size) ? 0 : 1;	// power of two
		// 	for (size_t r = z_size; r != 1; r /= 2) ++e;			// z_size <= 2^e
		// 	const int ln = e + 2;									// 64-bit => 16-bit
		// 	const size_t n = size_t(1) << ln;
		// 	if (ln > 32) { std::cout << z_size << ": size > 4G" << std::endl; break; }

		// 	// we must have n >= 16 * m_0 and m_0 = 2 * 4^e
		// 	const size_t m_0 = ((ln / 2) % 2 != 0) ? (size_t(1) << (ln / 2)) : (size_t(1) << (ln / 2 + 1));

		// 	// L1: 32 kB, L2: 256 kB / 2 MB
		// 	std::cout << z_size << ": " << ln << ", n = " << n << ", m_0 = " << m_0 << ", n / m_0 = " << n / m_0 <<
		// 		", m_0 size: " << m_0 * sizeof(Zp4) << std::endl;
		// }
		// return;

		// srand(time(nullptr));
		// for (size_t n = 7; n <= 20000000; n *= 2)
		// {
		// 	const size_t x_size = n, y_size = 2 * n / 3, z_size = x_size + y_size;
		// 	uint64_t * const x = new uint64_t[x_size];
		// 	uint64_t * const y = new uint64_t[y_size];
		// 	uint64_t * const z = new uint64_t[z_size];
		// 	uint64_t * const zr = new uint64_t[z_size];
		// 	Zp a = Zp(uint64_t(rand()));
		// 	for (size_t i = 0; i < x_size; ++i) { x[i] = a.get(); a *= Zp(55); }
		// 	for (size_t i = 0; i < y_size; ++i) { y[i] = a.get(); a *= Zp(55); }
		// 	for (size_t i = 0; i < z_size; ++i) z[i] = 0;
		// 	for (size_t i = 0; i < z_size; ++i) zr[i] = 0;

		// 	const size_t count = std::max(2000000 / n, size_t(1));

		// 	const auto tg = std::chrono::high_resolution_clock::now();
		// 	for (size_t i = 0; i < count; ++i) mpn_mul(mp_ptr(zr), mp_srcptr(x), mp_size_t(x_size), mp_srcptr(y), mp_size_t(y_size));
		// 	const double elapsed_time_g = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tg).count();

		// 	double elapsed_time_s = 0;
		// 	if (n < 4000)
		// 	{
		// 		const auto ts = std::chrono::high_resolution_clock::now();
		// 		for (size_t i = 0; i < count; ++i) smul(z, x, x_size, y, y_size);
		// 		elapsed_time_s = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - ts).count();
		// 		for (size_t i = 0; i < z_size; ++i) if (z[i] != zr[i]) std::cout << "Error!" << std::endl;
		// 	}

		// 	double elapsed_time_f = 0;
		// 	if (n > 128)
		// 	{
		// 		const auto tf = std::chrono::high_resolution_clock::now();
		// 		FastMul & fmul = FastMul::get_instance();
		// 		fmul.set_nthreads(4);
		// 		for (size_t i = 0; i < count; ++i)
		// 		{
		// 			fmul.init(z_size);
		// 			fmul.set_y(y, y_size);
		// 			fmul.mul_xy(x, x_size);
		// 			fmul.get_z(z, z_size);
		// 		}
		// 		elapsed_time_f = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - tf).count();
		// 		for (size_t i = 0; i < z_size; ++i) if (z[i] != zr[i]) std::cout << "Error!" << std::endl;
		// 	}

		// 	delete[] x;
		// 	delete[] y;
		// 	delete[] z;
		// 	delete[] zr;

		// 	std::cout << n << ", count: " << count << ", s: " << elapsed_time_s / elapsed_time_g << ", f: " << elapsed_time_f / elapsed_time_g
		// 			<< ", "<< elapsed_time_g << std::endl;
		// }
		// return;

		std::cout << header();
		if (argc < 2) std::cout << usage() << std::endl;

		const std::string arg1((argc > 1) ? argv[1] : "6912");

		int nthreads = 4;
		bool verbose = true;

		std::vector<std::string> args;
		for (int i = 1; i < argc; ++i) args.push_back(argv[i]);

		for (size_t i = 1, size = args.size(); i < size; ++i)
		{
			const std::string & arg = args[i];

			if (arg.substr(0, 2) == "-t")
			{
				const std::string ntstr = ((arg == "-t") && (i + 1 < size)) ? args[++i] : arg.substr(2);
				nthreads = std::atoi(ntstr.c_str());
			}
			if (arg == "-v") verbose = true;
		}

		CF & cf = CF::get_instance();
		cf.set_verbose(verbose);
		cf.set_nthreads(nthreads);

		guint N; N.from_string(arg1);
		if (!N.is_zero())
		{
			cf.solve(N);
		}
		else
		{
			// We must have N | 2^8 * 3^5 * 5^4 * 7^3 * 11^2 * 13^2 * 17^2 * 19^2 * 23 * ... * 199,
			// where the three dots represent the product of the primes between 23 and 199. See
			// P. Moree, H. J. J. te Riele, and J. Urbanowicz,
			// Divisibility properties of integers x, k satisfying 1^ + ... + (x − 1)^k = x^k,
			// Math. Comp. 63 (1994), 799-815.

			// N_max = 2^8 * 3^5 * 5^4 * 7^3
			static const std::vector<uint32_t> fN = { 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 5, 7, 7, 7 };

			for (size_t d = 0; d < fN.size(); ++d)
			{
				N = 1; for (size_t i = 0; i < d; ++i) N *= fN[i];
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
	}
};

int main(int argc, char * argv[])
{
	try
	{
		Application & app = Application::get_instance();
		app.run(argc, argv);
	}
	catch (const std::runtime_error & e)
	{
		std::cerr << std::endl << "Error: " << e.what() << "." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
