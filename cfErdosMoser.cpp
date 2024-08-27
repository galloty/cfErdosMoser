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
		ss << "cfErdosMoser 24.08.0 " << sysver << ssc.str() << std::endl;
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
		ss << " -t <n>: number of threads (1 or 4 threads, default: 4)," << std::endl;
		ss << " -v: verbose mode." << std::endl;
		return ss.str();
	}

public:
	void run(int argc, char * argv[])
	{
		std::cout << header();
		if (argc < 2) std::cout << usage() << std::endl;

		const std::string arg1((argc > 1) ? argv[1] : "6912");

		int nthreads = 4;
		bool verbose = false;

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
				if (!cf.solve(N)) return;
			}

			// Test all divisors of N_max
			std::set<uint64_t> Nset;
			for (uint64_t n2 = 1; n2 <= 256; n2 *= 2)
				for (uint64_t n3 = 1; n3 <= 243; n3 *= 3)
					for (uint64_t n5 = 1; n5 <= 625; n5 *= 5)
						for (uint64_t n7 = 1; n7 <= 343; n7 *= 7) Nset.insert(n2 * n3 * n5 * n7);
			Nset.erase(1);

			for (const uint64_t n : Nset) { N = n; if (!cf.solve(N)) return; }
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
	catch (...)
	{
		std::cerr << std::endl << "Error: unknown reason." << std::endl << std::flush;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
