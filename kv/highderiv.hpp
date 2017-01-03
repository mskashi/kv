/*
 * Copyright (c) 2015-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef HIGHDERIV_HPP
#define HIGHDERIV_HPP

#include <kv/psa.hpp>

namespace kv {

namespace ub = boost::numeric::ublas;

/*
 * calculate higher derivative of 1-variable function
 * using automatic differentiation
 *  - use interval as input if you want verified result.
 *  - if divfact == true then calculate f^{(i)}(x)/(i!)
 */

template <class F, class T>
ub::vector<T> highderiv(F f, const T& x, int n, bool divfact = false) {
	psa<T> p;
	int i;
	T fact;
	bool save_mode, save_uh, save_rh;

	save_mode = psa<T>::mode();
	save_uh = psa<T>::use_history();
	save_rh = psa<T>::record_history();
	psa<T>::mode() = 1;
	psa<T>::use_history() = false;
	psa<T>::record_history() = false;

	p.v.resize(2);
	p.v(0) = x;
	p.v(1) = 1;
	p = setorder(p, n);
	p = f(p);

	if (!divfact) {
		fact = 1.;
		for (i=1; i<p.v.size(); i++) {
			fact *= i;
			p.v(i) *= fact;
		}
	}

	psa<T>::mode() = save_mode;
	psa<T>::use_history() = save_uh;
	psa<T>::record_history() = save_rh;

	return p.v;
}

// make function object f^{n} from f
template <class F> struct HighDeriv {
	F f;
	int n;
	HighDeriv(F f, int n) : f(f), n(n) {}
	template <class T> T operator()(const T& x) {
		return kv::highderiv(f, x, n)[n];
	}
};

/*
 * fix except the n-th argument (partial application)
 */
template <class F, class TT>
struct Fix_Except_One_Arg {
	F f;
	int n;
	ub::vector<TT> x;
	int s;
	Fix_Except_One_Arg(F f, int n, const ub::vector<TT>& x): f(f), n(n), x(x) {
		s = x.size();
	}
	template <class T> T operator()(const T& in) {
		ub::vector<T> tmp(s);
		for (int i=0; i<s; i++) {
			if (i == n) {
				tmp(i) = in; // x(n) is discarded
			} else {
				tmp(i) = T(x(i));
			}
		}
		
		return f(tmp);
	}
};

// make [\partial^m f/\partial x_n^m] from f
template <class F> struct PartialDeriv {
	F f;
	int n, m;
	PartialDeriv(F f, int n, int m) : f(f), n(n), m(m) {}
	template <class T> T operator()(const ub::vector<T>& x) {
		Fix_Except_One_Arg<F,T> g(f, n, x);
		return kv::highderiv(g, x(n), m)[m];
	}
};

} // namespace kv;

#endif // HIGHDERIV_HPP
