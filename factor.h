/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

// Find the smallest factor of n using a sieve of Eratosthenes
class Factor
{
private:
	static const uint32_t sp_max = 1u << 19;
	static const size_t sieve_size = sp_max / 2;	// sieve using odd primes
	static const uint64_t n_max = sp_max * uint64_t(sp_max);

	std::vector<uint32_t> _prm;
	std::vector<size_t> _prm_ptr;
	uint32_t _factor[sieve_size];
	uint64_t _n, _j;

	void _step()
	{
		for (size_t k = 0; k < sieve_size; ++k) _factor[k] = 0;

		for (size_t i = 0, size = _prm.size(); i < size; ++i)
		{
			const uint32_t p = _prm[i];
			size_t k;
			for (k = _prm_ptr[i]; k < sieve_size; k += p) if (_factor[k] == 0) _factor[k] = p;
			_prm_ptr[i] = k - sieve_size;
		}
	}

	void _reset()
	{
		for (size_t i = 0, size = _prm.size(); i < size; ++i) _prm_ptr[i] = (_prm[i] >> 1) + size_t(_prm[i]);
		_step();
		_j = 0;
	}

public:
	void init() { _n = 1; _reset(); }

	Factor()
	{
		_prm.push_back(3); _prm.push_back(5); _prm.push_back(7);
		for (uint32_t k = 11; k < sp_max; k += 2)
		{
			const uint32_t s = uint32_t(std::sqrt(double(k))) + 1;
			uint32_t d;
			for (d = 3; d <= s; d += 2) if (k % d == 0) break;
			if (d > s) _prm.push_back(k);
		}
		_prm_ptr.resize(_prm.size());

		init();
	}

	virtual ~Factor() {}

	// smallest prime factor of n: faster if called in ascending order
	uint64_t smallest(const uint64_t n)
	{
		if (n % 2 == 0) return 2;
		if (n < _n) { std::cerr << "Warning: reset Factor." << std::endl; _reset(); }
		if (n >= n_max)
		{
			std::ostringstream ss; ss << "Factor max input is " << n_max;
			throw std::runtime_error(ss.str());
		}

		const uint64_t j = n / sp_max;
		_n = n;
		while (_j < j) { _step(); ++_j; }
		const uint64_t f = _factor[(n % sp_max) / 2];

		return (f == 0) ? n : f;
	}
};
