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
	Zp * _wi;
	Zp * _x;
	Zp * _y;

private:
	struct deleter { void operator()(const FastMul * const p) { delete p; } };

public:
	FastMul() : _n(0), _w(nullptr), _wi(nullptr), _x(nullptr), _y(nullptr) {}
	virtual ~FastMul() { free(); }
	
	static FastMul & get_instance()
	{
		static std::unique_ptr<FastMul, deleter> p_instance(new FastMul());
		return *p_instance;
	}

private:
	static size_t bitrev(const size_t i, const size_t n)
	{
		size_t r = 0;
		for (size_t k = n, j = i; k > 1; k /= 2, j /= 2) r = (2 * r) | (j % 2);
		return r;
	}

	void free()
	{
		Heap & heap = Heap::get_instance();
		if (_n != 0) { heap.free_fmul(_w, _n / 4); heap.free_fmul(_wi, _n / 4); heap.free_fmul(_x, _n); heap.free_fmul(_y, _n); }
		_n = 0; _w = _x = _y = nullptr;
	}

	void alloc(const size_t n)
	{
		free();
		Heap & heap = Heap::get_instance();
		_w = heap.alloc_fmul(n / 4); _wi = heap.alloc_fmul(n / 4); _x = heap.alloc_fmul(n); _y = heap.alloc_fmul(n);
		_n = n;
	}

	void check_size(const size_t n)
	{
		// Base = 2^16. We have (2^16 - 1)^2 * 2^32 < 2^64 - 2^32 + 1
		if (n > (size_t(1) << 32)) throw std::runtime_error("FastMul: size > 4G");

		// n_min = (16 + 16) * 4 = 128
		if (n > _n)
		{
			alloc(n);

			Zp * const w = _w;
			Zp * const wi = _wi;
			const Zp r = Zp::primroot_n(n), ri = r.invert();
			Zp t = Zp(1), ti = Zp(1);
			for (size_t i = 0; i < n / 2; ++i)
			{
				const size_t j = bitrev(i, n / 2);
				if (j < n / 4) { w[j] = t; wi[j] = ti; }
				t *= r; ti *= ri;
// std::cout << n << ", " << j << ": " << bitrev(j, n / 2) << ", " << w[j].get()<< ", " << (-wi[j]).get() << std::endl;
			}
// exit(0);
		}
	}

	void forward(const size_t n, const int ln, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const w = _w;

		for (size_t m = n / 4, s = n / 4 / m; m >= 2; m /= 4, s *= 4)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u2 = x[i + 2 * m], u1 = x[i + 1 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
				x[i + 0 * m] = v0 + v1; x[i + 1 * m] = v0 - v1; x[i + 2 * m] = v2 + v3; x[i + 3 * m] = v2 - v3;
			}

			for (size_t j = 1; j < s; ++j)
			{
				const Zp w_1 = w[j], w_2 = w[2 * j + 0], w_3 = w[2 * j + 1];	// w_1 = w_2 * w_2, w_3 = w_2.mul_i()

				for (size_t i = 0; i < m; ++i)
				{
					const size_t k = 4 * m * j + i;
					const Zp u0 = x[k + 0 * m], u2 = x[k + 2 * m] * w_1, u1 = x[k + 1 * m], u3 = x[k + 3 * m] * w_1;
					const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = (u1 + u3) * w_2, v3 = (u1 - u3) * w_3;
					x[k + 0 * m] = v0 + v1; x[k + 1 * m] = v0 - v1; x[k + 2 * m] = v2 + v3; x[k + 3 * m] = v2 - v3;
				}
			}
		}

		if (is_y)
		{
			if (ln % 2 == 0)
			{
				const Zp u0 = x[0], u2 = x[2], u1 = x[1], u3 = x[3];
				x[0] = u0 + u2; x[2] = u0 - u2; x[1] = u1 + u3; x[3] = u1 - u3;

				for (size_t j = 1; j < n / 4; ++j)
				{
					const Zp w_1 = w[j];

					const size_t k = 4 * j;
					const Zp u0 = x[k + 0], u2 = x[k + 2] * w_1, u1 = x[k + 1], u3 = x[k + 3] * w_1;
					x[k + 0] = u0 + u2; x[k + 2] = u0 - u2; x[k + 1] = u1 + u3; x[k + 3] = u1 - u3;
				}
			}
		}
	}

	void backward(const size_t n, const int ln, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const wi = _wi;

		for (size_t m = (ln % 2 == 0) ? 4 : 2, s = n / 4 / m; m <= n / 4; m *= 4, s /= 4)
		{
			for (size_t i = 0; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m], u2 = x[i + 2 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).mul_i();
				x[i + 0 * m] = v0 + v2; x[i + 2 * m] = v0 - v2; x[i + 1 * m] = v1 + v3; x[i + 3 * m] = v1 - v3;
			}

			for (size_t j = 1; j < s; ++j)
			{
				const Zp wi_1 = wi[j], wi_2 = wi[2 * j + 0], wi_3 = wi[2 * j + 1];

				for (size_t i = 0; i < m; ++i)
				{
					const size_t k = 4 * m * j + i;
					const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m], u2 = x[k + 2 * m], u3 = x[k + 3 * m];
					const Zp v0 = u0 + u1, v1 = (u0 - u1) * wi_2, v2 = u2 + u3, v3 = (u2 - u3) * wi_3;
					x[k + 0 * m] = v0 + v2; x[k + 2 * m] = (v0 - v2) * wi_1; x[k + 1 * m] = v1 + v3; x[k + 3 * m] = (v1 - v3) * wi_1;
				}
			}
		}
	}

	void mul(const size_t n, const int ln)
	{
		Zp * const x = _x;
		const Zp * const y = _y;
		const Zp * const w = _w;
		const Zp * const wi = _wi;

		if (ln % 2 == 0)
		{
			for (size_t j = 0; j < n / 4; ++j)
			{
				const Zp w_1 = w[j], wi_1 = wi[j];

				const size_t k = 4 * j;
				const Zp u0 = x[k + 0], u2 = x[k + 2] * w_1, u1 = x[k + 1], u3 = x[k + 3] * w_1;
				const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = u1 - u3;
				const Zp v0p = y[k + 0], v1p = y[k + 1], v2p = y[k + 2], v3p = y[k + 3];
				const Zp s0 = v0 * v0p + w_1 * v1 * v1p, s1 = v0 * v1p + v1 * v0p;
				const Zp s2 = v2 * v2p - w_1 * v3 * v3p, s3 = v2 * v3p + v3 * v2p;
				x[k + 0] = s0 + s2; x[k + 2] = (s0 - s2) * wi_1; x[k + 1] = s1 + s3; x[k + 3] = (s1 - s3) * wi_1;
			}
		}
		else
		{
			for (size_t j = 0; j < n / 4; ++j)
			{
				const Zp w_1 = w[j];

				const size_t k = 4 * j;
				const Zp u0 = x[k + 0], u1 = x[k + 1], u0p = y[k + 0], u1p = y[k + 1];
				x[k + 0] = u0 * u0p + w_1 * u1 * u1p; x[k + 1] = u0 * u1p + u1 * u0p;
				const Zp u2 = x[k + 2], u3 = x[k + 3], u2p = y[k + 2], u3p = y[k + 3];
				x[k + 2] = u2 * u2p - w_1 * u3 * u3p; x[k + 3] = u2 * u3p + u3 * u2p;
			}
		}
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

	static uint64_t get_d(const Zp * const x_i, const Zp & r, uint64_t & t)
	{
		uint64_t d_i = 0;
		for (size_t j = 0; j < 4; ++j)
		{
			t += Zp(x_i[j] * r).get();
			d_i |= uint64_t(uint16_t(t)) << (16 * j);
			t >>= 16;
		}
		return d_i;
	}

	void get(const size_t n, uint64_t * const d, const size_t d_size) const
	{
		const Zp * const x = _x;

		const Zp r = Zp::reciprocal(n / 2);

		uint64_t t0 = 0, t1 = 0;
		for (size_t i = 0, j = d_size / 2; i < d_size / 2; ++i, ++j)
		{
			d[i] = get_d(&x[4 * i], r, t0);
			d[j] = get_d(&x[4 * j], r, t1);
		}
		if (d_size % 2 == 1) d[d_size - 1] = get_d(&x[4 * (d_size - 1)], r, t1);
		for (size_t j = d_size / 2; j < d_size; ++j)
		{
			const __uint128_t t = d[j] + __uint128_t(t0);
			d[j] = uint64_t(t);
			t0 = uint64_t(t >> 64);
			if (t0 == 0) break;
		}
	}

public:
	void clear() { free(); }

	// x_size >= y_size >= 16, z_size = x_size + y_size
	void fmul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
	{
		const size_t z_size = x_size + y_size;
		int e = ((z_size & (~z_size + 1)) == z_size) ? 0 : 1;	// power of two
		for (size_t r = z_size; r != 1; r /= 2) ++e;			// z_size <= 2^e
		const int ln = e + 2; const size_t n = size_t(1) << ln;	// 64-bit => 16-bit

		check_size(n);

		set(n, x, x_size); set(n, y, y_size, true);

		forward(n, ln); forward(n, ln, true);
		mul(n, ln);
		backward(n, ln);

		get(n, z, x_size + y_size);
	}
};
