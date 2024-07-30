/*
Copyright 2024, Yves Gallot

cfErdosMoser is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#pragma once

#include <cstdint>
#include <stdexcept>
// #include <set>	// TODO remove

#include "mod64.h"
#include "heap.h"

class FastMul
{
private:
	// we must have n >= 16 * m_0
	static const size_t m_0 = 8 * 16;
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

		// static std::set<size_t> ns;
		// if (ns.find(n) == ns.end()) { std::cout << n << std::endl; ns.insert(n); }

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

	static void _forward(Zp * const x, const size_t m, const Zp & w_1, const Zp & w_2, const Zp & w_3)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const Zp u0 = x[i + 0 * m], u2 = x[i + 2 * m] * w_1, u1 = x[i + 1 * m], u3 = x[i + 3 * m] * w_1;
			const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = (u1 + u3) * w_2, v3 = (u1 - u3) * w_3;
			x[i + 0 * m] = v0 + v1; x[i + 1 * m] = v0 - v1; x[i + 2 * m] = v2 + v3; x[i + 3 * m] = v2 - v3;
		}
	}

	static void _backward(Zp * const x, const size_t m, const Zp & wi_1, const Zp & wi_2, const Zp & wi_3)
	{
		for (size_t i = 0; i < m; ++i)
		{
			const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m], u2 = x[i + 2 * m], u3 = x[i + 3 * m];
			const Zp v0 = u0 + u1, v1 = (u0 - u1) * wi_2, v2 = u2 + u3, v3 = (u2 - u3) * wi_3;
			x[i + 0 * m] = v0 + v2; x[i + 2 * m] = (v0 - v2) * wi_1; x[i + 1 * m] = v1 + v3; x[i + 3 * m] = (v1 - v3) * wi_1;
		}
	}

	static void _mul(Zp * const x, const Zp * const y, const Zp & w_1)
	{
			const Zp u0 = x[0], u1 = x[1], u0p = y[0], u1p = y[1];
			x[0] = u0 * u0p + w_1 * u1 * u1p; x[1] = u0 * u1p + u1 * u0p;
			const Zp u2 = x[2], u3 = x[3], u2p = y[2], u3p = y[3];
			x[2] = u2 * u2p - w_1 * u3 * u3p; x[3] = u2 * u3p + u3 * u2p;
	}

	void forward(const size_t n, const int ln, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const w = _w;

		if (ln % 2 == 0)
		{
			for (size_t i = 0, m = n / 8; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u4 = x[i + 4 * m], u2 = x[i + 2 * m], u6 = x[i + 6 * m];
				const Zp v0 = u0 + u4, v4 = u0 - u4, v2 = u2 + u6, v6 = Zp(u2 - u6).mul_i();
				const Zp u1 = x[i + 1 * m], u5 = x[i + 5 * m], u3 = x[i + 3 * m], u7 = x[i + 7 * m];
				const Zp v1 = u1 + u5, v5 = u1 - u5, v3 = u3 + u7, v7 = Zp(u3 - u7).mul_i();
				const Zp s0 = v0 + v2, s2 = v0 - v2, s1 = v1 + v3, s3 = Zp(v1 - v3).mul_i();
				x[i + 0 * m] = s0 + s1; x[i + 1 * m] = s0 - s1; x[i + 2 * m] = s2 + s3; x[i + 3 * m] = s2 - s3;
				const Zp s4 = v4 + v6, s6 = v4 - v6, s5 = Zp(v5 + v7).mul_sqrt_i(), s7 = Zp(v5 - v7).mul_i_sqrt_i();
				x[i + 4 * m] = s4 + s5; x[i + 5 * m] = s4 - s5; x[i + 6 * m] = s6 + s7; x[i + 7 * m] = s6 - s7;
			}
		}
		else
		{
			for (size_t i = 0, m = n / 4; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u2 = x[i + 2 * m], u1 = x[i + 1 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp(u1 - u3).mul_i();
				x[i + 0 * m] = v0 + v1; x[i + 1 * m] = v0 - v1; x[i + 2 * m] = v2 + v3; x[i + 3 * m] = v2 - v3;
			}
		}

		for (size_t m = (ln % 2 == 0) ? n / 32 : n / 16, s = n / 4 / m; m > m_0; m /= 4, s *= 4)
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
			for (size_t j_0 = 0, s_0 = n / 4 / m_0; j_0 < s_0; ++j_0)
			{
				for (size_t m = m_0, s = 1; m >= 2; m /= 4, s *= 4)
				{
					for (size_t j_s = 0; j_s < s; ++j_s)
					{
						const size_t j = j_0 * s + j_s;
						_forward(&x[4 * m * j], m, w[j], w[2 * j + 0], w[2 * j + 1]);
					}
				}
			}
		}
	}

	void backward(const size_t n, const int ln, const bool is_y = false)
	{
		Zp * const x = is_y ? _y : _x;
		const Zp * const wi = _wi;

		for (size_t m = 4 * m_0, s = n / 4 / m; m <= n / 16; m *= 4, s /= 4)
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

		if (ln % 2 == 0)
		{
			for (size_t i = 0, m = n / 8; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m], u2 = x[i + 2 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).mul_i();
				const Zp s0 = v0 + v2, s2 = v0 - v2, s1 = v1 + v3, s3 = v1 - v3;
				const Zp u4 = x[i + 4 * m], u5 = x[i + 5 * m], u6 = x[i + 6 * m], u7 = x[i + 7 * m];
				const Zp v4 = u4 + u5, v5 = Zp(u5 - u4).mul_i_sqrt_i(), v6 = u6 + u7, v7 = Zp(u7 - u6).mul_sqrt_i();
				const Zp s4 = v4 + v6, s6 = Zp(v6 - v4).mul_i(), s5 = v5 + v7, s7 = Zp(v7 - v5).mul_i();
				x[i + 0 * m] = s0 + s4; x[i + 4 * m] = s0 - s4; x[i + 2 * m] = s2 + s6; x[i + 6 * m] = s2 - s6;
				x[i + 1 * m] = s1 + s5; x[i + 5 * m] = s1 - s5; x[i + 3 * m] = s3 + s7; x[i + 7 * m] = s3 - s7;
			}
		}
		else
		{
			for (size_t i = 0, m = n / 4; i < m; ++i)
			{
				const Zp u0 = x[i + 0 * m], u1 = x[i + 1 * m], u2 = x[i + 2 * m], u3 = x[i + 3 * m];
				const Zp v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp(u3 - u2).mul_i();
				x[i + 0 * m] = v0 + v2; x[i + 2 * m] = v0 - v2; x[i + 1 * m] = v1 + v3; x[i + 3 * m] = v1 - v3;
			}
		}
	}

	void mul(const size_t n)
	{
		Zp * const x = _x;
		const Zp * const y = _y;
		const Zp * const w = _w;
		const Zp * const wi = _wi;

		// 4 * m_0 coefficients each step: 32 * m_0 bytes. 32 * 8 * 16 = 4 kB (L1: 32 kB, L2: 256 kB / 2 MB)
		for (size_t j_0 = 0, s_0 = n / 4 / m_0; j_0 < s_0; ++j_0)
		{
			for (size_t m = m_0, s = 1; m >= 2; m /= 4, s *= 4)
			{
				for (size_t j_s = 0; j_s < s; ++j_s)
				{
					const size_t j = j_0 * s + j_s;
					_forward(&x[4 * m * j], m, w[j], w[2 * j + 0], w[2 * j + 1]);	// w_1 = w_2 * w_2, w_3 = w_2.mul_i()
				}
			}

			for (size_t j_s = 0; j_s < m_0; ++j_s)
			{
				const size_t j = m_0 * j_0 + j_s;
				_mul(&x[4 * j], &y[4 * j], w[j]);
			}

			for (size_t m = 2, s = m_0 / m; m <= m_0; m *= 4, s /= 4)
			{
				for (size_t j_s = 0; j_s < s; ++j_s)
				{
					const size_t j = j_0 * s + j_s;
					_backward(&x[4 * m * j], m, wi[j], wi[2 * j + 0], wi[2 * j + 1]);
				}
			}
		}
	}

	static void set_d(Zp * const x_i, const uint64_t d_i)
	{
		uint64_t t = d_i;
		for (size_t j = 0; j < 4; ++j) { x_i[j] = Zp(uint16_t(t)); t >>= 16; }
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

	void set(const size_t n, const uint64_t * const d, const size_t d_size, const bool is_y = false)
	{
		Zp * const x = is_y ? _y :_x;

		for (size_t i = 0; i < d_size; ++i) set_d(&x[4 * i], d[i]);
		for (size_t i = 4 * d_size; i < n; ++i) x[i] = Zp(0);
	}

	void get(const size_t n, uint64_t * const d, const size_t d_size) const
	{
		const Zp * const x = _x;

		const Zp r = Zp::reciprocal(n / 2);
		uint64_t t = 0; for (size_t i = 0; i < d_size; ++i) d[i] = get_d(&x[4 * i], r, t);
	}

public:
	void clear() { free(); }

	// x_size >= y_size >= 129, z_size = x_size + y_size
	void fmul(uint64_t * const z, const uint64_t * const x, const size_t x_size, const uint64_t * const y, const size_t y_size)
	{
		const size_t z_size = x_size + y_size;
		int e = ((z_size & (~z_size + 1)) == z_size) ? 0 : 1;	// power of two
		for (size_t r = z_size; r != 1; r /= 2) ++e;			// z_size <= 2^e
		const int ln = e + 2; const size_t n = size_t(1) << ln;	// 64-bit => 16-bit

		// n >= (x_size + y_size) * 4 = 2 * 129 * 4 = 1032 => n_min = 2048
		check_size(n);

		set(n, x, x_size); set(n, y, y_size, true);

		forward(n, ln); forward(n, ln, true);
		mul(n);
		backward(n, ln);

		get(n, z, x_size + y_size);
	}
};
