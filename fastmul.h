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
	int _ln;
	size_t _m_0, _ni_4;
	size_t _wsize, _xsize;
	Zp * _w;
	Zp * _wi;
	Zp4 * _x;
	Zp4 * _y;

private:
	struct deleter { void operator()(const FastMul * const p) { delete p; } };

public:
	FastMul() : _ln(0), _m_0(0), _ni_4(0), _wsize(0), _xsize(0), _w(nullptr), _wi(nullptr), _x(nullptr), _y(nullptr) {}
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

	static const size_t _index_gap = 64 / sizeof(Zp4);	// Cache line size is 64 bytes
	static constexpr size_t index(const size_t k, const size_t m_0) { const size_t j = k / m_0, i = k % m_0; return j * (m_0 + _index_gap) + i; }

	void free()
	{
		Heap & heap = Heap::get_instance();
		if (_wsize != 0)
		{
			heap.free_fmul(_w, _wsize * sizeof(Zp));
			heap.free_fmul(_wi, _wsize * sizeof(Zp));
			heap.free_fmul(_x, _xsize * sizeof(Zp4));
			heap.free_fmul(_y, _xsize * sizeof(Zp4));
		}
		_wsize = _xsize = 0; _w = _wi = nullptr; _x = _y = nullptr;
	}

	void alloc(const size_t n)
	{
		free();
		Heap & heap = Heap::get_instance();
		_wsize = n / 4; _xsize = _ni_4;
		_w = static_cast<Zp *>(heap.alloc_fmul(_wsize * sizeof(Zp)));
		_wi = static_cast<Zp *>(heap.alloc_fmul(_wsize * sizeof(Zp)));
		_x = static_cast<Zp4 *>(heap.alloc_fmul(_xsize * sizeof(Zp4)));
		_y = static_cast<Zp4 *>(heap.alloc_fmul(_xsize * sizeof(Zp4)));
	}

	finline static void _vforward4(Zp4 * const x, const size_t m, const Zp & w_1, const Zp & w_2, const Zp & w_3)
	{
		const Zp4 u0 = x[0 * m], u2 = x[2 * m] * w_1, u1 = x[1 * m], u3 = x[3 * m] * w_1;
		const Zp4 v0 = u0 + u2, v2 = u0 - u2, v1 = (u1 + u3) * w_2, v3 = (u1 - u3) * w_3;
		x[0 * m] = v0 + v1; x[1 * m] = v0 - v1; x[2 * m] = v2 + v3; x[3 * m] = v2 - v3;
	}

	finline static void _vbackward4(Zp4 * const x, const size_t m, const Zp & wi_1, const Zp & wi_2, const Zp & wi_3)
	{
		const Zp4 u0 = x[0 * m], u1 = x[1 * m], u2 = x[2 * m], u3 = x[3 * m];
		const Zp4 v0 = u0 + u1, v1 = (u0 - u1) * wi_2, v2 = u2 + u3, v3 = (u2 - u3) * wi_3;
		x[0 * m] = v0 + v2; x[2 * m] = (v0 - v2) * wi_1; x[1 * m] = v1 + v3; x[3 * m] = (v1 - v3) * wi_1;
	}

	finline static void _vforward4v2(Zp4 * const x, const Zp2 & w_1, const Zp2 & w_2, const Zp2 & w_3)
	{
		Zp4::shuffle2in(x);
		const Zp4 u0 = x[0], u2 = x[2] * w_1, u1 = x[1], u3 = x[3] * w_1;
		const Zp4 v0 = u0 + u2, v2 = u0 - u2, v1 = (u1 + u3) * w_2, v3 = (u1 - u3) * w_3;
		x[0] = v0 + v1; x[1] = v0 - v1; x[2] = v2 + v3; x[3] = v2 - v3;
		Zp4::shuffle2out(x);
	}

	finline static void _vbackward4v2(Zp4 * const x, const Zp2 & wi_1, const Zp2 & wi_2, const Zp2 & wi_3)
	{
		Zp4::shuffle2in(x);
		const Zp4 u0 = x[0], u1 = x[1], u2 = x[2], u3 = x[3];
		const Zp4 v0 = u0 + u1, v1 = (u0 - u1) * wi_2, v2 = u2 + u3, v3 = (u2 - u3) * wi_3;
		x[0] = v0 + v2; x[2] = (v0 - v2) * wi_1; x[1] = v1 + v3; x[3] = (v1 - v3) * wi_1;
		Zp4::shuffle2out(x);
	}

	finline static void _vforward4_0(Zp4 * const x, const size_t m)
	{
		const Zp4 u0 = x[0 * m], u2 = x[2 * m], u1 = x[1 * m], u3 = x[3 * m];
		const Zp4 v0 = u0 + u2, v2 = u0 - u2, v1 = u1 + u3, v3 = Zp4(u1 - u3).mul_i();
		x[0 * m] = v0 + v1; x[1 * m] = v0 - v1; x[2 * m] = v2 + v3; x[3 * m] = v2 - v3;
	}

	finline static void _vbackward4_0(Zp4 * const x, const size_t m)
	{
		const Zp4 u0 = x[0 * m], u1 = x[1 * m], u2 = x[2 * m], u3 = x[3 * m];
		const Zp4 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp4(u3 - u2).mul_i();
		x[0 * m] = v0 + v2; x[2 * m] = v0 - v2; x[1 * m] = v1 + v3; x[3 * m] = v1 - v3;
	}

	finline static void _vforward8_0(Zp4 * const x, const size_t m)
	{
		const Zp4 u0 = x[0 * m], u4 = x[4 * m], u2 = x[2 * m], u6 = x[6 * m];
		const Zp4 v0 = u0 + u4, v4 = u0 - u4, v2 = u2 + u6, v6 = Zp4(u2 - u6).mul_i();
		const Zp4 u1 = x[1 * m], u5 = x[5 * m], u3 = x[3 * m], u7 = x[7 * m];
		const Zp4 v1 = u1 + u5, v5 = u1 - u5, v3 = u3 + u7, v7 = Zp4(u3 - u7).mul_i();
		const Zp4 s0 = v0 + v2, s2 = v0 - v2, s1 = v1 + v3, s3 = Zp4(v1 - v3).mul_i();
		x[0 * m] = s0 + s1; x[1 * m] = s0 - s1; x[2 * m] = s2 + s3; x[3 * m] = s2 - s3;
		const Zp4 s4 = v4 + v6, s6 = v4 - v6, s5 = Zp4(v5 + v7).mul_sqrt_i(), s7 = Zp4(v5 - v7).mul_i_sqrt_i();
		x[4 * m] = s4 + s5; x[5 * m] = s4 - s5; x[6 * m] = s6 + s7; x[7 * m] = s6 - s7;
	}

	finline static void _vbackward8_0(Zp4 * const x, const size_t m)
	{
		const Zp4 u0 = x[0 * m], u1 = x[1 * m], u2 = x[2 * m], u3 = x[3 * m];
		const Zp4 v0 = u0 + u1, v1 = u0 - u1, v2 = u2 + u3, v3 = Zp4(u3 - u2).mul_i();
		const Zp4 s0 = v0 + v2, s2 = v0 - v2, s1 = v1 + v3, s3 = v1 - v3;
		const Zp4 u4 = x[4 * m], u5 = x[5 * m], u6 = x[6 * m], u7 = x[7 * m];
		const Zp4 v4 = u4 + u5, v5 = Zp4(u5 - u4).mul_i_sqrt_i(), v6 = u6 + u7, v7 = Zp4(u7 - u6).mul_sqrt_i();
		const Zp4 s4 = v4 + v6, s6 = Zp4(v6 - v4).mul_i(), s5 = v5 + v7, s7 = Zp4(v7 - v5).mul_i();
		x[0 * m] = s0 + s4; x[4 * m] = s0 - s4; x[2 * m] = s2 + s6; x[6 * m] = s2 - s6;
		x[1 * m] = s1 + s5; x[5 * m] = s1 - s5; x[3 * m] = s3 + s7; x[7 * m] = s3 - s7;
	}

	finline static void _mul4(Zp4 & x, const Zp4 & y, const Zp & w_1)
	{
		// const Zp u0 = x[0], u1 = x[1], u2 = x[2], u3 = x[3], u0p = y[0], u1p = y[1], u2p = y[2], u3p = y[3];
		// const Zp v0 = u0 * u0p + w_1 * u1 * u1p, v1 = u0 * u1p + u1 * u0p;
		// const Zp v2 = u2 * u2p - w_1 * u3 * u3p, v3 = u2 * u3p + u3 * u2p;
		// x = Zp4(v0, v1, v2, v3);

		const Zp4 xy = x * y;				// u0 * u0p, u1 * u1p, u2 * u2p, u3 * u3p
		const Zp4 xyp = x * y.permute();	// u0 * u1p, u1 * u0p, u2 * u3p, u3 * u2p
		const Zp4 e = xy.even(xyp);			// u0 * u0p, u0 * u1p, u2 * u2p, u2 * u3p
		const Zp4 o = xy.odd(xyp);			// u1 * u1p, u1 * u0p, u3 * u3p, u3 * u2p
		const Zp4 op = o.mul_w(w_1);		// w_1 * u1 * u1p, u1 * u0p, -w_1 * u3 * u3p, u3 * u2p
		x = e + op;
	}

	void forward(const bool is_y = false)
	{
		const int ln = _ln; const size_t n = size_t(1) << ln;
		const size_t m_0 = _m_0, mi_0 = m_0 + _index_gap, ni_4 = _ni_4;
		Zp4 * const x = is_y ? _y : _x;
		const Zp * const w = _w;

		// n / (2 * m_0) Zp4 each step
		for (size_t i_t = 0; i_t < m_0 / 2; ++i_t)
		{
			if (ln % 2 == 0)
			{
				for (size_t i_m = 2 * i_t, mi = ni_4 / 8; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vforward8_0(&x[i_m + i], mi);
				}
			}
			else
			{
				for (size_t i_m = 2 * i_t, mi = ni_4 / 4; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vforward4_0(&x[i_m + i], mi);
				}
			}

			for (size_t s = (ln % 2 == 0) ? 8 : 4, mi = ni_4 / (4 * s); mi >= mi_0; mi /= 4, s *= 4)
			{
				for (size_t i_m = 2 * i_t; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vforward4_0(&x[i_m + i], mi);
				}

				for (size_t j = 1; j < s; ++j)
				{
					const Zp w_1 = w[j], w_2 = w[2 * j + 0], w_3 = w[2 * j + 1];	// w_1 = w_2 * w_2, w_3 = w_2.mul_i()

					for (size_t i_m = 2 * i_t; i_m < mi; i_m += mi_0)
					{
						for (size_t i = 0; i < 2; ++i) _vforward4(&x[4 * mi * j + i_m + i], mi, w_1, w_2, w_3);
					}
				}
			}
		}

		if (is_y)
		{
			for (size_t j_t = 0, s_0 = n / 4 / m_0; j_t < s_0; ++j_t)
			{
				Zp4 * const xt = &x[(m_0 + _index_gap) * j_t];

				for (size_t m = m_0 / 4, s = 1; m >= 2; m /= 4, s *= 4)
				{
					for (size_t j_s = 0; j_s < s; ++j_s)
					{
						const size_t j = j_t * s + j_s;
						const Zp w_1 = w[j], w_2 = w[2 * j + 0], w_3 = w[2 * j + 1];
						for (size_t i = 0; i < m; ++i) _vforward4(&xt[4 * m * j_s + i], m, w_1, w_2, w_3);
					}
				}

				for (size_t j_s = 0, s = m_0 / 4; j_s < s; ++j_s)
				{
					const size_t j = j_t * s + j_s;
					const Zp2 w_1 = Zp2(w[2 * j + 0], w[2 * j + 1]), w_2 = Zp2(w[4 * j + 0], w[4 * j + 2]), w_3 = Zp2(w[4 * j + 1], w[4 * j + 3]);
					_vforward4v2(&xt[4 * j_s], w_1, w_2, w_3);
				}
			}
		}
	}

	void backward()
	{
		const int ln = _ln; const size_t n = size_t(1) << ln;
		const size_t m_0 = _m_0, mi_0 = m_0 + _index_gap, ni_4 = _ni_4;
		Zp4 * const x = _x;
		const Zp * const wi = _wi;

		// n / (2 * m_0) Zp4 each step
		for (size_t i_t = 0; i_t < 4 * m_0 / 8; ++i_t)
		{
			for (size_t s = n / 16 / m_0, mi = mi_0; mi <= ni_4 / 16; mi *= 4, s /= 4)
			{
				for (size_t i_m = 2 * i_t; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vbackward4_0(&x[i_m + i], mi);
				}

				for (size_t j = 1; j < s; ++j)
				{
					const Zp wi_1 = wi[j], wi_2 = wi[2 * j + 0], wi_3 = wi[2 * j + 1];

					for (size_t i_m = 2 * i_t; i_m < mi; i_m += mi_0)
					{
						for (size_t i = 0; i < 2; ++i) _vbackward4(&x[4 * mi * j + i_m + i], mi, wi_1, wi_2, wi_3);
					}
				}
			}

			if (ln % 2 == 0)
			{
				for (size_t i_m = 2 * i_t, mi = ni_4 / 8; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vbackward8_0(&x[i_m + i], mi);
				}
			}
			else
			{
				for (size_t i_m = 2 * i_t, mi = ni_4 / 4; i_m < mi; i_m += mi_0)
				{
					for (size_t i = 0; i < 2; ++i) _vbackward4_0(&x[i_m + i], mi);
				}
			}
		}
	}

	void mult()
	{
		const size_t n = size_t(1) << _ln;
		const size_t m_0 = _m_0;
		Zp4 * const x = _x;
		const Zp4 * const y = _y;
		const Zp * const w = _w;
		const Zp * const wi = _wi;

		// m_0 consecutive Zp4 each step
		for (size_t j_t = 0, s_0 = n / 4 / m_0; j_t < s_0; ++j_t)
		{
			const size_t k_t = (m_0 + _index_gap) * j_t;
			Zp4 * const xt = &x[k_t];
			const Zp4 * const yt = &y[k_t];

			for (size_t m = m_0 / 4, s = 1; m >= 2; m /= 4, s *= 4)
			{
				for (size_t j_s = 0; j_s < s; ++j_s)
				{
					const size_t j = j_t * s + j_s;
					const Zp w_1 = w[j], w_2 = w[2 * j + 0], w_3 = w[2 * j + 1];
					for (size_t i = 0; i < m; ++i) _vforward4(&xt[4 * m * j_s + i], m, w_1, w_2, w_3);	// w_1 = w_2 * w_2, w_3 = w_2.mul_i()
				}
			}

			for (size_t j_s = 0, s = m_0 / 4; j_s < s; ++j_s)
			{
				const size_t j = j_t * s + j_s;
				const Zp2 w_1 = Zp2(w[2 * j + 0], w[2 * j + 1]), w_2 = Zp2(w[4 * j + 0], w[4 * j + 2]), w_3 = Zp2(w[4 * j + 1], w[4 * j + 3]);
				_vforward4v2(&xt[4 * j_s], w_1, w_2, w_3);
			}

			for (size_t j_s = 0; j_s < m_0; ++j_s)
			{
				const size_t j = m_0 * j_t + j_s;
				_mul4(xt[j_s], yt[j_s], w[j]);
			}

			for (size_t j_s = 0, s = m_0 / 4; j_s < s; ++j_s)
			{
				const size_t j = j_t * s + j_s;
				const Zp2 wi_1 = Zp2(wi[2 * j + 0], wi[2 * j + 1]), wi_2 = Zp2(wi[4 * j + 0], wi[4 * j + 2]), wi_3 = Zp2(wi[4 * j + 1], wi[4 * j + 3]);
				_vbackward4v2(&xt[4 * j_s], wi_1, wi_2, wi_3);
			}

			for (size_t m = 2, s = m_0 / (4 * m); m <= m_0 / 4; m *= 4, s /= 4)
			{
				for (size_t j_s = 0; j_s < s; ++j_s)
				{
					const size_t j = j_t * s + j_s;
					const Zp wi_1 = wi[j], wi_2 = wi[2 * j + 0], wi_3 = wi[2 * j + 1];
					for (size_t i = 0; i < m; ++i) _vbackward4(&xt[4 * m * j_s + i], m, wi_1, wi_2, wi_3);
				}
			}
		}
	}

	void set(const uint64_t * const d, const size_t d_size, const bool is_y = false)
	{
		const size_t m_0 = _m_0, ni_4 = _ni_4;
		Zp4 * const x = is_y ? _y :_x;

		size_t k = 0, ki = 0, i = 0;
		while (k < d_size)
		{
			const uint64_t d_k = d[k];
			x[ki] = Zp4(Zp(uint16_t(d_k)), Zp(uint16_t(d_k >> 16)), Zp(uint16_t(d_k >> 32)), Zp(uint16_t(d_k >> 48)));

			++i; ++k; ++ki;
			if (i == m_0) { i = 0; ki += _index_gap; }
		}
		while (ki < ni_4) { x[ki] = Zp4(Zp(0)); ++ki; }
	}

	void get(uint64_t * const d, const size_t d_size) const
	{
		const size_t m_0 = _m_0;
		const Zp4 * const x = _x;

		const Zp r = Zp::reciprocal(_ln - 1);
		__uint128_t t = 0;
		size_t k = 0, ki = 0, i = 0;
		while (k < d_size)
		{
			const Zp4 xr = x[ki] * r;
			uint64_t d_k = 0;
			for (size_t j = 0; j < 4; ++j)
			{
				t += xr[j].get();
				d_k |= uint64_t(uint16_t(t)) << (16 * j);
				t >>= 16;
			}
			d[k] = d_k;

			++i; ++k; ++ki;
			if (i == m_0) { i = 0; ki += _index_gap; }
		}
	}

public:
	void clear() { free(); }

	// x_size >= y_size >= 129, z_size = x_size + y_size
	void init(const size_t z_size)
	{
		int e = ((z_size & (~z_size + 1)) == z_size) ? 0 : 1;	// power of two
		for (size_t r = z_size; r != 1; r /= 2) ++e;			// z_size <= 2^e
		const int ln = e + 2;									// 64-bit => 16-bit
		_ln = ln; const size_t n = size_t(1) << ln;
		// we must have n >= 16 * m_0 and m_0 = 2 * 4^h
		_m_0 = ((ln / 2) % 2 != 0) ? (size_t(1) << (ln / 2)) : (size_t(1) << (ln / 2 + 1));
		_ni_4 = index(n / 4, _m_0);

		// Base = 2^16. We have (2^16 - 1)^2 * 2^32 < 2^64 - 2^32 + 1
		if (ln > 32) throw std::runtime_error("FastMul: size > 4G");

		if (n / 4 > _wsize)
		{
			alloc(n);

			Zp * const w = _w;
			Zp * const wi = _wi;
			const Zp r = Zp::primroot_n(ln), ri = r.invert();
			Zp t = Zp(1), ti = Zp(1);
			for (size_t i = 0; i < n / 2; ++i)
			{
				const size_t j = bitrev(i, n / 2);
				if (j < n / 4) { w[j] = t; wi[j] = ti; }
				t *= r; ti *= ri;
			}
		}
	}

	void set_x(const uint64_t * const d, const size_t d_size) { set(d, d_size); }
	void set_y(const uint64_t * const d, const size_t d_size) { set(d, d_size, true); }
	void get_x(uint64_t * const d, const size_t d_size) { get(d, d_size); }
	void mul() { forward(); forward(true); mult(); backward(); }
};
