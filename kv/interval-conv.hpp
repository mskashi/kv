/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef INTERVAL_CONV_HPP
#define INTERVAL_CONV_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

namespace kv {

template <int N>
void mpfrtodouble(const mpfr<N>& x, double& y, int rnd = 0)
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
void ddtompfr(const dd& x, mpfr<N>& y, int rnd = 0)
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

template <int N>
void mpfrtodd(const mpfr<N>& x, dd& y, int rnd = 0)
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

	mpfrtodouble(x, dtmp1, 0);

	// mtmp = x - dtmp1;
	// theoretically no rounding error.
	// use rounded subtraction just to be safe.
	mpfr_sub_d(mtmp.a, x.a, dtmp1, mode);

	dtmp2 = mpfr_get_d(mtmp.a, mode);

	kv::dd::twosum(dtmp1, dtmp2, y.a1, y.a2);
}

template <int N, int M>
void mpfrtompfr(const mpfr<N>& x, mpfr<M>& y, int rnd = 0)
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

void ddtodouble(const dd& x, double& y, int rnd = 0)
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

void iddtoidouble(const interval<dd>& x, interval<double>& y)
{
	ddtodouble(x.lower(), y.lower(), -1);
	ddtodouble(x.upper(), y.upper(), 1);
}

void idoubletoidd(const interval<double>& x, interval<dd>& y)
{
	y.lower() = x.lower();
	y.upper() = x.upper();
}

template <int N>
void impfrtoidouble(const interval< mpfr<N> >& x, interval<double>& y)
{
	mpfrtodouble(x.lower(), y.lower(), -1);
	mpfrtodouble(x.upper(), y.upper(), 1);
}

template <int N>
void idoubletoimpfr(const interval<double>& x, interval< mpfr<N> >& y)
{
	y.lower() = x.lower();
	y.upper() = x.upper();
}

template <int N>
void impfrtoidd(const interval< mpfr<N> >& x, interval<dd>& y)
{
	mpfrtodd(x.lower(), y.lower(), -1);
	mpfrtodd(x.upper(), y.upper(), 1);
}

template <int N>
void iddtoimpfr(const interval<dd>& x, interval< mpfr<N> >& y)
{
	ddtompfr(x.lower(), y.lower(), -1);
	ddtompfr(x.upper(), y.upper(), 1);
}

template <int N, int M>
void impfrtoimpfr(const interval< mpfr<N> >& x, interval< mpfr<M> >& y)
{
	mpfrtompfr(x.lower(), y.lower(), -1);
	mpfrtompfr(x.upper(), y.upper(), 1);
}

};

#endif // INTERVAL_CONV_HPP
