/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>

#include "mod64.h"
#include "heap.h"

class FastMul
{
private:
	size_t _n;
	Zp * _w;
	Zp * _x;
	Zp * _y;

private:
	struct deleter { void operator()(const FastMul * const p) { delete p; } };

public:
	FastMul() : _n(0), _w(nullptr), _x(nullptr), _y(nullptr) {}
	virtual ~FastMul() {}

	static FastMul & get_instance()
	{
		static std::unique_ptr<FastMul, deleter> p_instance(new FastMul());
		return *p_instance;
	}

private:
	void check_size(const size_t n)
	{
		// Base = 2^16. We have (2^16 - 1)^2 * 2^32 < 2^64 - 2^32 + 1
		if (n > (size_t(1) << 32)) throw std::runtime_error("FastMul: size > 4G");

		if (n > _n)
		{
			Heap & heap = Heap::get_instance();
			if (_n != 0) { heap.free_fmul(_w, _n); heap.free_fmul(_x, _n); heap.free_fmul(_y, _n); }

			_n = n;
			_w = heap.allocate_fmul(n); _x = heap.allocate_fmul(n); _y = heap.allocate_fmul(n);

			Zp * const w = _w;
			for (size_t m = 1; m <= n / 2; m *= 2)
			{
				const Zp r = Zp::primroot_n(2 * m);
				Zp t = Zp(1);
				for (size_t j = 0; j < m; ++j) { w[m + j] = t; t *= r; }
			}
		}
	}

	void forward(const size_t n, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const w = _w;

		for (size_t s = 1, m = n / 2; m >= 1; s *= 2, m /= 2)
		{
			for (size_t j = 0; j < s; ++j)
			{
				const size_t k = 2 * j * m;
				const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
				x[k + 0 * m] = u0 + u1; x[k + 1 * m] = u0 - u1;

				for (size_t i = 1; i < m; ++i)
				{
					const Zp w_si = w[m + i];
					const size_t k = 2 * j * m + i;
					const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
					x[k + 0 * m] = u0 + u1;
					x[k + 1 * m] = (u0 - u1) * w_si;
				}
			}
		}
	}

	void backward(const size_t n)
	{
		Zp * const x = _x;
		const Zp * const w = _w;

		for (size_t s = n / 2, m = 1; s >= 1; s /= 2, m *= 2)
		{
			for (size_t j = 0; j < s; ++j)
			{
				const size_t k = 2 * j * m;
				const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m];
				x[k + 0 * m] = u0 + u1; x[k + 1 * m] = u0 - u1;

				for (size_t i = 1; i < m; ++i)
				{
					const Zp w_si_inv = -w[m + m - i];
					const size_t k = 2 * j * m + i;
					const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m] * w_si_inv;
					x[k + 0 * m] = u0 + u1;
					x[k + 1 * m] = u0 - u1;
				}
			}
		}

		const Zp r = Zp::reciprocal(n);
		for (size_t k = 0; k < n; ++k) x[k] *= r;
	}

	void mul(const size_t n)
	{
		Zp * const x = _x;
		const Zp * const y = _y;
		for (size_t k = 0; k < n; ++k) x[k] *= y[k];
	}

	void sqr(const size_t n)
	{
		Zp * const x = _x;
		for (size_t k = 0; k < n; ++k) x[k] *= x[k];
	}

	void set(const size_t n, const uint64_t * const d, const size_t d_size, const bool is_y = false)
	{
		Zp * const x = is_y ? _y :_x;

		for (size_t i = 0; i < d_size; ++i)
		{
			uint64_t d_i = d[i];
			for (size_t j = 0; j < 4; ++j) { x[4 * i + j] = Zp(uint16_t(d_i)); d_i >>= 16; }
		}
		for (size_t i = 4 * d_size; i < n; ++i) x[i] = Zp(0);
	}

	void get(uint64_t * const d, const size_t d_size) const
	{
		const Zp * const x = _x;

		uint64_t t = 0;
		for (size_t i = 0; i < d_size; ++i)
		{
			uint64_t d_i[4];
			for (size_t j = 0; j < 4; ++j) { t += x[4 * i + j].get(); d_i[j] = uint16_t(t); t >>= 16; }
			d[i] = d_i[0] | (d_i[1] << 16) | (d_i[2] << 32) | (d_i[3] << 48);
		}
	}

public:
	void clear()
	{
		if (_n != 0) { Heap & heap = Heap::get_instance(); heap.free_fmul(_w, _n); heap.free_fmul(_x, _n); heap.free_fmul(_y, _n); }
		_n = 0; _w = _x = _y = nullptr;
	}

	// x_size >= y_size >= 16, z_size = x_size + y_size
	void fmul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
	{
		size_t e = ((x_size & (~x_size + 1)) == x_size) ? 0 : 1;	// power of two
		for (size_t r = x_size; r != 1; r /= 2) ++e;				// x_size <= 2^e

		const size_t n = size_t(1) << (e + 1 + 2);					// +1: twice the size, +2: 64-bit => 16-bit

		check_size(n);

		set(n, x, x_size); set(n, y, y_size, true);

		forward(n); forward(n, true);
		mul(n);
		backward(n);

		get(z, x_size + y_size);
	}
};
