#ifndef STROBOMAP_HPP
#define STROBOMAP_HPP

#include <exception>
#include "ode-nv.hpp"
#include "ode-autodif-nv.hpp"
#include "ode-maffine.hpp"
#ifdef USE_MAFFINE2
#include "ode-maffine2.hpp"
#endif

namespace ub = boost::numeric::ublas;


namespace kv {


// 微分方程式の右辺関数fの関数オブジェクトから、
// ストロボマップの関数オブジェクトを生成する。
// 生成された関数オブジェクトは、引数として、
//   vector<T>,
//   vector< autodif<T> >,
//   vector< interval<T> >,
//   vector< affine<T> >,
//   vector< autodif< interval<T> > >
// を受け付ける。それぞれ、内部で
//   odelong_nv (通常版),
//   odelong_nv (autodif版),
//   odelong_maffine (通常版),
//   odelong_maffine (affine版),
//   odelong_maffine (autodif版)
// を呼び出している。(-DUSE_MAFFINE2ならmaffine2)

template <class F, class TT> class StroboMap {
	public:
	F f;
	TT start, end;
	int order;
	int iter_max;
	int verbose;
	int ep_reduce;
	int ep_reduce_limit;

	StroboMap(F f_v, TT start_v, TT end_v, int order_v, int iter_max_v = 2, int verbose_v = 0, int ep_reduce_v = 0, int ep_reduce_limit_v = 0)
	: f(f_v), start(start_v), end(end_v), order(order_v), iter_max(iter_max_v), verbose(verbose_v), ep_reduce(ep_reduce_v), ep_reduce_limit(ep_reduce_limit_v) {
	}

	// mid_ifnecessary<T1,T2>(x)は、T2からT1に変換可能ならxそのまま、
	// 不可能ならmid(x)を返す。

	#include <boost/utility/enable_if.hpp>

	template <class T1, class T2> T1 inline static mid_ifnecessary(T2& x, typename boost::enable_if_c< convertible<T2, T1>::value >::type* =0) {
		return T1(x);
	}

	template <class T1, class T2> T1 inline static mid_ifnecessary(T2& x, typename boost::enable_if_c< ! convertible<T2, T1>::value >::type* =0) {
		return T1(mid(x));
	}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> result;
		T start2, end2;

		start2 = mid_ifnecessary<T,TT>(start);
		end2 = mid_ifnecessary<T,TT>(end);

		result = x;

		odelong_nv(f, result, start2, end2, order, verbose);

		return result;
	}

	template <class T> ub::vector< autodif<T> > operator() (const ub::vector< autodif<T> >& x){
		ub::vector< autodif<T> > result;
		T start2, end2;

		start2 = mid_ifnecessary<T,TT>(start);
		end2 = mid_ifnecessary<T,TT>(end);

		result = x;

		odelong_nv(f, result, start2, end2, order, verbose);

		return result;
	}

	template <class T> ub::vector< interval<T> > operator() (const ub::vector< interval<T> >& x){
		ub::vector< interval<T> > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		#ifdef USE_MAFFINE2
		r = odelong_maffine2(f, result, (interval<T>)start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
		#else
		r = odelong_maffine(f, result, (interval<T>)start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
		#endif

		if (r != 2) {
			throw std::range_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}

	template <class T> ub::vector< affine<T> > operator() (const ub::vector< affine<T> >& x){
		ub::vector< affine<T> > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		#ifdef USE_MAFFINE2
		r = odelong_maffine2(f, result, (interval<T>)start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
		#else
		r = odelong_maffine(f, result, (interval<T>)start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
		#endif

		if (r != 2) {
			throw std::range_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}

	template <class T> ub::vector< autodif< interval<T> > > operator() (const ub::vector< autodif< interval<T> > >& x){
		ub::vector< autodif< interval<T> > > result;
		interval<T> end2;
		int r;

		result = x;
		end2 = end;

		r = odelong_maffine(f, result, (interval<T>)start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
		if (r != 2) {
			throw std::range_error("StroboMap(): cannot calculate validated solution.");
		}

		return result;
	}
};

// おまけ:
// 関数fの関数オブジェクトから関数x-f(x)の関数オブジェクトを生成

template <class F> class FixedPoint {
	public:
	F f;
	FixedPoint(F f_v) : f(f_v) {}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		return x - f(x);
	}
};


// 2変数関数のstrobomapから、2点境界値問題をShooting法で解くための
// 非線形方程式を生成
//   x(start_index) = start_value, x(end_index) = end_value
//   のような境界条件を指定する。
//   x(start_indexでないindex)を未知数として変動させ、
//   e(end_index) - end_value を0にしようとする。

template <class F, class TV> class Shooting_TPBVP {
	public:
	F f;
	TV start_value, end_value;
	int start_index, end_index;
	int variable_index;

	Shooting_TPBVP(F f_v, TV start_value_v, TV end_value_v, int start_index_v, int end_index_v) : f(f_v) , start_value(start_value_v), end_value(end_value_v), start_index(start_index_v), end_index(end_index_v) {
		variable_index = 1 - start_index;
	}

	// mid_ifnecessary<T1,T2>(x)は、T2からT1に変換可能ならxそのまま、
	// 不可能ならmid(x)を返す。

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
