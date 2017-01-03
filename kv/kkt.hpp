/*
 * Copyright (c) 2015-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef KKT_HPP
#define KKT_HPP

/*
 * make KKT equation using automatic differentiation
 */

#include <kv/interval.hpp>
#include <kv/autodif.hpp>

namespace kv {

namespace ub = boost::numeric::ublas;


// alpha_plus/alpha_minus function for non-interval value

template <class T>
T alpha_plus(const T& x) {
	using std::pow;
	if (x <= 0.) return 0.;
	else return pow(x, 3);
}

template <class T>
T alpha_minus(const T& x) {
	using std::pow;
	if (x >= 0.) return 0.;
	else return -pow(x, 3);
}

template <class T>
T alpha_plus_d(const T& x) {
	using std::pow;
	if (x <= 0.) return 0.;
	else return 3 * pow(x, 2);
}

template <class T>
T alpha_minus_d(const T& x) {
	using std::pow;
	if (x >= 0.) return 0.;
	else return -3 * pow(x, 2);
}


// alpha_plus/alpha_minus function for interval value

template <class T>
interval<T> alpha_plus(const interval<T>& x) {
	interval<T> r;

	if (x.lower() <= 0) {
		r.lower() = 0.;
	} else {
		r.lower() = pow(interval<T>(x.lower()), 3.).lower();
	}

	if (x.upper() <= 0) {
		r.upper() = 0.;
	} else {
		r.upper() = pow(interval<T>(x.upper()), 3.).upper();
	}

	return r;
}

template <class T>
interval<T> alpha_plus_d(const interval<T>& x) {
	interval<T> r;

	if (x.lower() <= 0) {
		r.lower() = 0.;
	} else {
		r.lower() = (3. * pow(interval<T>(x.lower()), 2.)).lower();
	}

	if (x.upper() <= 0) {
		r.upper() = 0.;
	} else {
		r.upper() = (3. * pow(interval<T>(x.upper()), 2.)).upper();
	}

	return r;
}

template <class T>
interval<T> alpha_minus(const interval<T>& x) {
	interval<T> r;

	if (x.lower() >= 0) {
		r.upper() = 0.;
	} else {
		r.upper() = (-pow(interval<T>(x.lower()), 3)).upper();
	}

	if (x.upper() >= 0) {
		r.lower() = 0.;
	} else {
		r.lower() = (-pow(interval<T>(x.upper()), 3)).lower();
	}

	return r;
}

template <class T>
interval<T> alpha_minus_d(const interval<T>& x) {
	interval<T> r;

	if (x.lower() >= 0) {
		r.lower() = 0.;
	} else {
		r.lower() = (-3. * pow(interval<T>(x.lower()), 2)).lower();
	}

	if (x.upper() >= 0) {
		r.upper() = 0.;
	} else {
		r.upper() = (-3. * pow(interval<T>(x.upper()), 2)).upper();
	}

	return r;
}

// alpha_plus/alpha_minus function for autodif

template <class T>
autodif<T> alpha_plus(const autodif<T>& x) {
	autodif<T> r;

	r.v = alpha_plus(x.v);
	r.d = alpha_plus_d(x.v) * x.d;

	return r;
}

template <class T>
autodif<T> alpha_minus(const autodif<T>& x) {
	autodif<T> r;

	r.v = alpha_minus(x.v);
	r.d = alpha_minus_d(x.v) * x.d;

	return r;
}

// KKT equation maker using objective function, inequalities, equalities

template <class F1, class F2, class F3>
struct KKT_equation {
	F1 f;
	F2 g;
	F3 h;

	KKT_equation(F1 f, F2 g, F3 h) : f(f), g(g), h(h) {}

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		int sf = f.size();
		int sg, sh;
		ub::vector<T> xx;
		ub::vector<T> fv;
		ub::vector<T> gv;
		ub::matrix<T> gm;
		ub::vector<T> hv;
		ub::matrix<T> hm;
		ub::vector<T> b;
		ub::vector<T> mu;
		ub::vector<T> r;
		int i;
		T dummy;

		xx.resize(sf);
		for (i=0; i<sf; i++) xx(i) = x(i);

		autodif<T>::split(f(autodif<T>::init(xx)), dummy, fv);

		autodif<T>::split(g(autodif<T>::init(xx)), gv, gm);
		sg = gv.size();
		if (sg != 0) {
			b.resize(sg);
			for (i=0; i<sg; i++) {
				b(i) = alpha_plus(x(sf + i));
				gv(i) += alpha_minus(x(sf + i));
			}
			fv += prod(trans(gm), b);
		}

		autodif<T>::split(h(autodif<T>::init(xx)), hv, hm);
		sh = hv.size();
		if (sh != 0) {
			mu.resize(sh);
			for (i=0; i<sh; i++) mu(i) = x(sf + sg + i);
			fv += prod(trans(hm), mu);
		}

		r.resize(sf + sg + sh);
		for (i=0; i<sf; i++) r(i) = fv(i);
		for (i=0; i<sg; i++) r(sf + i) = gv(i);
		for (i=0; i<sh; i++) r(sf + sg + i) = hv(i);

		return r;
	}
};

}; // namespace kv

#endif // KKT_HPP
