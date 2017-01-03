/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
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

} // namespace kv;

#endif // HIGHDERIV_HPP
