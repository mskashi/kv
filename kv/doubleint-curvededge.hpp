/*
 * Copyright (c) 2021 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DOUBLEINT_CURVEDEDGE
#define DOUBLEINT_CURVEDEDGE

#include <iostream>
#include "double-newtoncotes.hpp"

namespace kv {

template <class F1, class F2, class F3>
struct Convert_Curved_Edge {
	F1 f;
	F2 el;
	F3 eh;
	Convert_Curved_Edge(F1 f, F2 el, F3 eh) : f(f), el(el), eh(eh) {}
	template <class T> T operator()(T x, T y) {
		T tmp, tmp2;
		tmp = el(x);
		tmp2 = eh(x) - tmp;
		return f(x, tmp + y * tmp2) * tmp2;
	}
};

template <class F>
struct Swap_xy {
	F f;
	Swap_xy(F f) : f(f) {}
	template <class T> T operator()(T x, T y) {
		return f(y, x);
	}
};

/*
 *  \int_start^end ( \int_el(x)^eh(x) f(x,y) dy ) dx
 */

template <class T, class F1, class F2, class F3>
interval<T>
doubleint_curvededge(F1 f, F2 el, F3 eh, interval<T> start, interval<T> end, int n, int m1, int m2, int errorterm_div = 1)
{
	Convert_Curved_Edge<F1, F2, F3> g(f, el, eh);
	return double_newtoncotes(g, start, end, interval<T>(0), interval<T>(1), n, m1, m2, errorterm_div);
}


/*
 *  \int_start^end ( \int_el(y)^eh(y) f(x,y) dx ) dy
 */

template <class T, class F1, class F2, class F3>
interval<T>
doubleint_curvededge_s(F1 f, F2 el, F3 eh, interval<T> start, interval<T> end, int n, int m1, int m2, int errorterm_div = 1)
{
	Swap_xy<F1> g(f);
	return doubleint_curvededge(g, el, eh, start, end, n, m2, m1, errorterm_div);
}

}; // namespace kv

#endif // DOUBLEINT_CURVEDEDGE
