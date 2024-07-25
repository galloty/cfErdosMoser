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

		// n_min = (16 + 16) * 4 = 128
		if (n > _n)
		{
			Heap & heap = Heap::get_instance();
			if (_n != 0) { heap.free_fmul(_w, _n); heap.free_fmul(_x, _n); heap.free_fmul(_y, _n); }

			_n = n;
			_w = heap.allocate_fmul(n); _x = heap.allocate_fmul(n); _y = heap.allocate_fmul(n);

			Zp * const w = _w;
			for (size_t m = 2; m <= n / 4; m *= 2)
			{
				const Zp r = Zp::primroot_n(4 * m), r_inv = r.invert();
				Zp t = r, t_inv = r_inv;
				for (size_t i = 1; i < m; ++i)
				{
					w[2 * m + i] = t; t *= r;
					w[2 * m + m + i] = t_inv; t_inv *= r_inv;
				}
			}
		}
	}

	void forward(const size_t n, const int ln, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const w = _w;

		for (size_t m = n / 4, s = n / 4 / m; m >= 2; m /= 4, s *= 4)
		{
			for (size_t j = 0; j < s; ++j)
			{
				const size_t k = 4 * m * j;
				const Zp u0 = x[k + 0 * m], u2 = x[k + 2 * m], u1 = x[k + 1 * m], u3 = x[k + 3 * m];
				const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
				x[k + 0 * m] = v0 + v1; x[k + 1 * m] = v0 - v1; x[k + 2 * m] = v2 + v3; x[k + 3 * m] = v2 - v3;

				for (size_t i = 1; i < m; ++i)
				{
					const Zp w_2 = w[2 * m + i], w_1 = w_2 * w_2, w_3 = w_1 * w_2;
					const size_t k = 4 * m * j + i;
					const Zp u0 = x[k + 0 * m], u2 = x[k + 2 * m], u1 = x[k + 1 * m], u3 = x[k + 3 * m];
					const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
					x[k + 0 * m] = v0 + v1; x[k + 1 * m] = (v0 - v1) * w_1; x[k + 2 * m] = (v2 + v3) * w_2; x[k + 3 * m] = (v2 - v3) * w_3;
				}
			}
		}

		if (is_y)
		{
			if (ln % 2 == 0)
			{
				for (size_t k = 0; k < n; k += 4)
				{
					const Zp u0 = x[k + 0], u2 = x[k + 2], u1 = x[k + 1], u3 = x[k + 3];
					const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
					x[k + 0] = v0 + v1; x[k + 1] = v0 - v1; x[k + 2] = v2 + v3; x[k + 3] = v2 - v3;
				}
			}
			else
			{
				for (size_t k = 0; k < n; k += 4)
				{
					const Zp u0 = x[k + 0], u1 = x[k + 1], u2 = x[k + 2], u3 = x[k + 3];
					x[k + 0] = u0 + u1; x[k + 1] = u0 - u1; x[k + 2] = u2 + u3; x[k + 3] = u2 - u3;
				}
			}
		}
	}

	void backward(const size_t n, const int ln)
	{
		Zp * const x = _x;
		const Zp * const w = _w;

		for (size_t m = (ln % 2 == 0) ? 4 : 2, s = n / 4 / m; m <= n / 4; m *= 4, s /= 4)
		{
			for (size_t j = 0; j < s; ++j)
			{
				const size_t k = 4 * m * j;
				const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m], u2 = x[k + 2 * m], u3 = x[k + 3 * m];
				const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).mul_i();
				x[k + 0 * m] = v0 + v2; x[k + 2 * m] = v0 - v2; x[k + 1 * m] = v1 + v3; x[k + 3 * m] = v1 - v3;

				for (size_t i = 1; i < m; ++i)
				{
					const Zp w_2i = w[2 * m + m + i], w_1i = w_2i * w_2i, w_3i = w_1i * w_2i;
					const size_t k = 4 * m * j + i;
					const Zp u0 = x[k + 0 * m], u1 = x[k + 1 * m] * w_1i, u2 = x[k + 2 * m] * w_2i, u3 = x[k + 3 * m] * w_3i;
					const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).mul_i();
					x[k + 0 * m] = v0 + v2; x[k + 2 * m] = v0 - v2; x[k + 1 * m] = v1 + v3; x[k + 3 * m] = v1 - v3;
				}
			}
		}
	}

	void mul(const size_t n, const int ln)
	{
		Zp * const x = _x;
		const Zp * const y = _y;

		if (ln % 2 == 0)
		{
			for (size_t k = 0; k < n; k += 4)
			{
				const Zp u0 = x[k + 0], u2 = x[k + 2], u1 = x[k + 1], u3 = x[k + 3];
				const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
				const Zp s0 = (v0 + v1) * y[k + 0], s1 = (v0 - v1) * y[k + 1];
				const Zp s2 = (v2 + v3) * y[k + 2], s3 = (v2 - v3) * y[k + 3];
				const Zp t0 = s0 + s1, t1 = s0 - s1, t2 = s2 + s3, t3 = Zp(s3 - s2).mul_i();
				x[k + 0] = t0 + t2; x[k + 2] = t0 - t2; x[k + 1] = t1 + t3; x[k + 3] = t1 - t3;
			}
		}
		else
		{
			for (size_t k = 0; k < n; k += 4)
			{
				const Zp u0 = x[k + 0], u1 = x[k + 1], u2 = x[k + 2], u3 = x[k + 3];
				const Zp v0 = (u0 + u1) * y[k + 0], v1 = (u0 - u1) * y[k + 1];
				const Zp v2 = (u2 + u3) * y[k + 2], v3 = (u2 - u3) * y[k + 3];
				x[k + 0] = v0 + v1; x[k + 1] = v0 - v1; x[k + 2] = v2 + v3; x[k + 3] = v2 - v3;
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

	void get(const size_t n, uint64_t * const d, const size_t d_size) const
	{
		const Zp * const x = _x;

		const Zp r = Zp::reciprocal(n);

		uint64_t t = 0;
		for (size_t i = 0; i < d_size; ++i)
		{
			uint64_t d_i[4];
			for (size_t j = 0; j < 4; ++j) { const Zp x_ij = x[4 * i + j] * r; t += x_ij.get(); d_i[j] = uint16_t(t); t >>= 16; }
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
