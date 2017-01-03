/*
 * Copyright (c) 2013-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef STROBOMAP_HPP
#define STROBOMAP_HPP

#include <stdexcept>
#include <kv/ode-nv.hpp>
#include <kv/ode-autodif-nv.hpp>
#include <kv/ode-maffine.hpp>
#ifdef USE_MAFFINE2
#include <kv/ode-maffine2.hpp>
#endif
#include <kv/ode-param.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;


// Generate function object of strobomap
// using function object of r.h.s of differential equation.
// Generated function object can receive
//   vector<T>,
//   vector< autodif<T> >,
//   vector< interval<T> >,
//   vector< affine<T> >,
//   vector< autodif< interval<T> > >
// as argument. In each case,
//   odelong_nv in ode-nv.hpp,
//   odelong_nv in ode-autodif-nv.hpp,
//   odelong_maffine in ode-maffine.hpp (interval version)
//   odelong_maffine in ode-maffine.hpp (affine version),
//   odelong_maffine in ode-maffine.hpp (autodif version)
// is called inside.
// (If -DUSE_MAFFINE2 then ode-maffine2.hpp is used instead.)

template <class F, class T> class StroboMap {
	public:
	F f;
	interval<T> start, end;
	ode_param<T> p;

	StroboMap(F f, interval<T> start, interval<T> end, ode_param<T> p = ode_param<T>())
	: f(f), start(start), end(end), p(p) {}

	ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> result;

		result = x;

		odelong_nv(f, result, mid(start), mid(end), p);

		return result;
	}

	ub::vector< autodif<T> > operator() (const ub::vector< autodif<T> >& x){
		ub::vector< autodif<T> > result;

		result = x;

		odelong_nv(f, result, mid(start), mid(end), p);

		return result;
	}

	ub::vector< interval<T> > operator() (const ub::vector< interval<T> >& x){
		ub::vector< interval<T> > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		#ifdef USE_MAFFINE2
		r = odelong_maffine2(f, result, start, end2, p);
		#else
		r = odelong_maffine(f, result, start, end2, p);
		#endif

		if (r != 2) {
			throw std::domain_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}

	ub::vector< affine<T> > operator() (const ub::vector< affine<T> >& x){
		ub::vector< affine<T> > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		#ifdef USE_MAFFINE2
		r = odelong_maffine2(f, result, start, end2, p);
		#else
		r = odelong_maffine(f, result, start, end2, p);
		#endif

		if (r != 2) {
			throw std::domain_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}

	ub::vector< autodif< interval<T> > > operator() (const ub::vector< autodif< interval<T> > >& x){
		ub::vector< autodif< interval<T> > > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		r = odelong_maffine(f, result, start, end2, p);
		if (r != 2) {
			throw std::domain_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}
};

// Generate function object of "x-f(x)" from function object of f.

template <class F> class FixedPoint {
	public:
	F f;
	FixedPoint(F f) : f(f) {}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		return x - f(x);
	}
};


// Generate function object for shooting method of two point
// boundary value problem from 2-variable strobomap.
//   boundary condition:
//     x(start_index) = start_value, x(end_index) = end_value
//   unknown variable:
//     x(index != start_index)
// Solve ODE from initial value (start_value, unknown) and try to
// satisfy x(end_index) = end_value by changing unknown.

template <class F, class TV> class Shooting_TPBVP {
	public:
	F f;
	TV start_value, end_value;
	int start_index, end_index;
	int variable_index;

	Shooting_TPBVP(F f, TV start_value, TV end_value, int start_index, int end_index) : f(f) , start_value(start_value), end_value(end_value), start_index(start_index), end_index(end_index) {
		variable_index = 1 - start_index;
	}

	// mid_ifnecessary<T1,T2>(x) returns x if T2 is convertible to T1,
	// and returns mid(x) if impossoble.

	#include <boost/utility/enable_if.hpp>

	template <class T1, class T2> T1 inline static mid_ifnecessary(T2& x, typename boost::enable_if_c< convertible<T2, T1>::value >::type* =0) {
		return T1(x);
	}

	template <class T1, class T2> T1 inline static mid_ifnecessary(T2& x, typename boost::enable_if_c< ! convertible<T2, T1>::value >::type* =0) {
		return T1(mid(x));
	}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& s) {
		ub::vector<T> x, y, r;
		x.resize(2);
		r.resize(1);
		x(variable_index) = s(0);
		x(start_index) = mid_ifnecessary<T,TV>(start_value);
		y = f(x);
		r(0) = y(end_index) - mid_ifnecessary<T,TV>(end_value);
		return r;
	}
};

} // namespace kv

#endif // STROBOMAP_HPP
