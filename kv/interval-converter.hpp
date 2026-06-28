/*
 * Copyright (c) 2022-2026 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef INTERVAL_CONVERTER_HPP
#define INTERVAL_CONVERTER_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

namespace kv {

// dd to double

inline void rounded_converter(const dd& x, double& y, int rnd = 0)
{
	if (rnd == 1) {
		rop<double>::begin();
		y = rop<double>::add_up(x.a1, x.a2);
		rop<double>::end();
	} else if (rnd == -1) {
		rop<double>::begin();
		y = rop<double>::add_down(x.a1, x.a2);
		rop<double>::end();
	} else {
		y = x.a1 + x.a2;
	}
}

// double to dd

inline void rounded_converter(const double& x, dd& y, int rnd = 0)
{
	y = x;
}

// mpfr to (double, dd)

template <int N>
void rounded_converter(const mpfr<N>& x, double& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	y = mpfr_get_d(x.a, mode);
}

template <int N>
void rounded_converter(const mpfr<N>& x, dd& y, int rnd = 0)
{
	mpfr<N> mtmp;
	double dtmp1, dtmp2;

	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	rounded_converter(x, dtmp1, 0);

	// mtmp = x - dtmp1;
	// theoretically no rounding error.
	// use rounded subtraction just to be safe.
	mpfr_sub_d(mtmp.a, x.a, dtmp1, mode);

	dtmp2 = mpfr_get_d(mtmp.a, mode);

	dd::twosum(dtmp1, dtmp2, y.a1, y.a2);
}

// (double, dd) to mpfr

template <int N>
void rounded_converter(const double& x, mpfr<N>& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	mpfr_set_d(y.a, x, mode);
}

template <int N>
void rounded_converter(const dd& x, mpfr<N>& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	// y = x.a;
	mpfr_set_d(y.a, x.a1, mode);
	mpfr_add_d(y.a, y.a, x.a2, mode);
}


// mpfr<N> to mpfr<M>

template <int N, int M>
void rounded_converter(const mpfr<N>& x, mpfr<M>& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	mpfr_set(y.a, x.a, mode);
}

}; // namespace kv;


// converters for the types using fp80

#include <kv/fp80.hpp>

#ifdef KV_HAVE_FP80

#include <kv/rfp80.hpp>
#include <kv/ddx.hpp>
#include <kv/rddx.hpp>

#include <kv/hwround.hpp>
#include <limits>
#include <cmath>

namespace kv {

// fp80 to (double, dd, mpfr)

inline void rounded_converter(const fp80& x, double& y, int rnd = 0)
{
	volatile fp80 x1 = x;
	volatile double r;

	if (rnd == 1) {
		hwround::roundup();
	} else if (rnd == -1) {
		hwround::rounddown();
	}

	r = x1;

	if (rnd == 1 || rnd == -1) {
		hwround::roundnear();
	}

	y = r;
}

inline void rounded_converter(const fp80& x, dd& y, int rnd = 0)
{
	fp80 tmp;
	double z1, z2;

	rounded_converter(x, z1, 0);
	if (z1 == std::numeric_limits<double>::infinity()) {
		if (rnd == -1) {
			y = (std::numeric_limits<dd>::max)();
		} else {
			y = std::numeric_limits<dd>::infinity();
		}
		return;
	}
	if (z1 == -std::numeric_limits<double>::infinity()) {
		if (rnd == 1) {
			y = -(std::numeric_limits<dd>::max)();
		} else {
			y = -std::numeric_limits<dd>::infinity();
		}
		return;
	}

	// error free?
	tmp = x - z1;

	rounded_converter(tmp, z2, rnd);

	dd::twosum(z1, z2, y.a1, y.a2);
	if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
		y.a2 = 0;
	}
}

template <int N>
void rounded_converter(const fp80& x, mpfr<N>& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	mpfr_set_ld(y.a, x, mode);
}

// (double, dd, mpfr) to fp80

inline void rounded_converter(const double& x, fp80& y, int rnd = 0)
{
	y = x;
}

inline void rounded_converter(const dd& x, fp80& y, int rnd = 0)
{
	fp80 z1 = x.a1, z2 = x.a2;

	if (rnd == 1) {
		rop<fp80>::begin();
		y = rop<fp80>::add_up(z1, z2);
		rop<fp80>::end();
	} else if (rnd == -1) {
		rop<fp80>::begin();
		y = rop<fp80>::add_down(z1, z2);
		rop<fp80>::end();
	} else {
		y = z1 + z2;
	}
}

template <int N>
void rounded_converter(const mpfr<N>& x, fp80& y, int rnd = 0)
{
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	y = mpfr_get_ld(x.a, mode);
}

// ddx to (double, dd, mpfr, fp80)

inline void rounded_converter(const ddx& x, double& y, int rnd = 0)
{
	volatile fp80 x1 = x.a1, x2 = x.a2;
	volatile double r;

	if (rnd == 1) {
		hwround::roundup();
	} else if (rnd == -1) {
		hwround::rounddown();
	}

	r = x1 + x2;

	if (rnd == 1 || rnd == -1) {
		hwround::roundnear();
	}

	y = r;
}

inline void rounded_converter(const ddx& x, dd& y, int rnd = 0)
{
	fp80 tmp;
	double z1, z2;

	rounded_converter(x.a1, z1, 0);
	if (z1 == std::numeric_limits<double>::infinity()) {
		if (rnd == -1) {
			y = (std::numeric_limits<dd>::max)();
		} else {
			y = std::numeric_limits<dd>::infinity();
		}
		return;
	}
	if (z1 == -std::numeric_limits<double>::infinity()) {
		if (rnd == 1) {
			y = -(std::numeric_limits<dd>::max)();
		} else {
			y = -std::numeric_limits<dd>::infinity();
		}
		return;
	}

	// error free?
	tmp = x.a1 - z1;

	if (rnd == 1) {
		rop<fp80>::begin();
		tmp = rop<fp80>::add_up(tmp, x.a2);
		rop<fp80>::end();
	} else if (rnd == -1) {
		rop<fp80>::begin();
		tmp = rop<fp80>::add_down(tmp, x.a2);
		rop<fp80>::end();
	} else {
		tmp += x.a2;
	}

	rounded_converter(tmp, z2, rnd);

	dd::twosum(z1, z2, y.a1, y.a2);
	if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
		y.a2 = 0;
	}
}

template <int N>
void rounded_converter(const ddx& x, mpfr<N>& y, int rnd = 0)
{
	mpfr<N> tmp;
	mp_rnd_t mode;

	if (rnd == 1) {
		mode = MPFR_RNDU;
	} else if (rnd == -1) {
		mode = MPFR_RNDD;
	} else {
		mode = MPFR_RNDN;
	}

	mpfr_set_ld(y.a, x.a1, mode);
	mpfr_set_ld(tmp.a, x.a2, mode);
	mpfr_add(y.a, y.a, tmp.a, mode);
}

inline void rounded_converter(const ddx& x, fp80& y, int rnd = 0)
{
	if (rnd == 1) {
		rop<fp80>::begin();
		y = rop<fp80>::add_up(x.a1, x.a2);
		rop<fp80>::end();
	} else if (rnd == -1) {
		rop<fp80>::begin();
		y = rop<fp80>::add_down(x.a1, x.a2);
		rop<fp80>::end();
	} else {
		y = x.a1 + x.a2;
	}
}

// (double, dd, mpfr, fp80) to ddx

inline void rounded_converter(const double& x, ddx& y, int rnd = 0)
{
	y = x;
}

inline void rounded_converter(const dd& x, ddx& y, int rnd = 0)
{
	fp80 dtmp1 = x.a1, dtmp2 = x.a2;

	ddx::twosum(dtmp1, dtmp2, y.a1, y.a2);
}

template <int N>
void rounded_converter(const mpfr<N>& x, ddx& y, int rnd = 0)
{
	mpfr<N> mtmp1, mtmp2;
	fp80 dtmp1, dtmp2;

	rounded_converter(x, dtmp1, 0);
	mtmp1 = x;
	rounded_converter(dtmp1, mtmp2);
	mtmp1 -= mtmp2;
	rounded_converter(mtmp1, dtmp2, rnd);

	ddx::twosum(dtmp1, dtmp2, y.a1, y.a2);
}

inline void rounded_converter(const fp80& x, ddx& y, int rnd = 0)
{
	y = x;
}

}; // namespace kv;

#endif // KV_HAVE_FP80

namespace kv {

// interval converter
// called by operator= in interval.hpp

template <class T1, class T2>
void interval_converter(const interval<T1>& x, interval<T2>& y)
{
	rounded_converter(x.lower(), y.lower(), -1);
	rounded_converter(x.upper(), y.upper(), 1);
}

};

#endif // INTERVAL_CONVERTER_HPP
