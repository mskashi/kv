/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef AUTODIF_HPP
#define AUTODIF_HPP

// Automatic Differentiation by bottom up algorithm

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

#include <kv/convert.hpp>


namespace ub = boost::numeric::ublas;


namespace kv {

	template <class T> class autodif;
	template <class C, class T> struct convertible<C, autodif<T> > {
		static const bool value = convertible<C, T>::value || boost::is_same<C, autodif<T> >::value;
	};
	template <class C, class T> struct acceptable_n<C, autodif<T> > {
		static const bool value = convertible<C, T>::value;
	};


template <class T> class autodif {
	public:
	T v;
	ub::vector<T> d;

	typedef T base_type;

	autodif() {
		v = 0.;
		d.resize(0);
	}

	template <class C> explicit autodif(const C& x, typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value >::type* =0) {
		v = x;
		d.resize(0);
	}

	template <class C> typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif& >::type operator=(const C& x) {
		v = x;
		d.resize(0);
		return *this;
	}

	friend autodif operator+(const autodif& a, const autodif& b) {
		autodif r;

		r.v = a.v + b.v;

		if (a.d.size() == 0) {
			r.d = b.d;
		} else if (b.d.size() == 0) {
			r.d = a.d;
		} else {
			r.d = a.d + b.d;
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator+(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v + b;
		r.d = a.d;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator+(const C& a, const autodif& b) {
		autodif r;

		r.v = a + b.v;
		r.d = b.d;

		return r;
	}

	friend autodif& operator+=(autodif& a, const autodif& b) {
		a = a + b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif& >::type operator+=(autodif& a, const C& b) {
		a.v += b;
		return a;
	}

	friend autodif operator-(const autodif& a, const autodif& b) {
		autodif r;

		r.v = a.v - b.v;

		if (a.d.size() == 0) {
			r.d = - b.d;
		} else if (b.d.size() == 0) {
			r.d = a.d;
		} else {
			r.d = a.d - b.d;
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator-(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v - b;
		r.d = a.d;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator-(const C& a, const autodif& b) {
		autodif r;

		r.v = a - b.v;
		r.d = - b.d;

		return r;
	}

	friend autodif& operator-=(autodif& a, const autodif& b) {
		a = a - b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif& >::type operator-=(autodif& a, const C& b) {
		a.v -= b;
		return a;
	}

	friend autodif operator-(const autodif& a) {
		autodif r;

		r.v = - a.v;
		r.d = - a.d;

		return r;
	}

	friend autodif operator*(const autodif& a, const autodif& b) {
		autodif r;

		r.v = a.v * b.v;

		if (a.d.size() == 0) {
			r.d = a.v * b.d;
		} else if (b.d.size() == 0) {
			r.d = b.v * a.d;
		} else {
			r.d = b.v * a.d + a.v * b.d;
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator*(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v * b;
		// r.d = b * a.d;
		r.d = T(b) * a.d; // assist for VC++

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator*(const C& a, const autodif& b) {
		autodif r;

		r.v = a * b.v;
		// r.d = a * b.d;
		r.d = T(a) * b.d; // assist for VC++

		return r;
	}

	friend autodif& operator*=(autodif& a, const autodif& b) {
		a = a * b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif& >::type operator*=(autodif& a, const C& b) {
		a.v *= b;
		// a.r *= b;
		a.r *= T(b); // assist for VC++
		return a;
	}

	friend autodif operator/(const autodif& a, const autodif& b) {
		autodif r;

		r.v = a.v / b.v;

		if (a.d.size() == 0) {
			r.d = b.d * (-a.v/(b.v*b.v));
		} else if (b.d.size() == 0) {
			r.d = a.d / b.v;
		} else {
			r.d = a.d / b.v + b.d * (-a.v/(b.v*b.v));
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator/(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v / b;
		// r.d = a.d / b;
		r.d = a.d / T(b); // assist for VC++

		return r;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif >::type operator/(const C& a, const autodif& b) {
		autodif r;

		r.v = a / b.v;
		r.d = b.d * (-a/(b.v*b.v));

		return r;
	}

	friend autodif& operator/=(autodif& a, const autodif& b) {
		a = a / b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< kv::acceptable_n<C, autodif>::value, autodif& >::type operator/=(autodif& a, const C& b) {
		a.v /= b;
		// a.d /= b;
		a.d /= T(b); // assist for VC++
		return a;
	}

	friend std::ostream& operator<<(std::ostream& s, const autodif& x) {
		int i;
		int n = x.d.size();
		s << x.v;
		s << '<';
		for (i=0; i<n; i++) {
			s << x.d(i);
			if (i != n-1) {
				s << ',';
			}
		}
		s << '>';
		return s;
	}


	friend autodif pow(const autodif& x, int y) {
		autodif r;

		using std::pow;
		r.v = pow(x.v, y);
		r.d = (y * pow(x.v, y - 1)) * x.d;
		return r;
	}

	friend autodif pow(const autodif& x, const autodif& y) {
		return exp(y * log(x));
	}


	friend autodif exp (const autodif& x) {
		autodif r;

		using std::exp;
		r.v = exp(x.v);
		// r.d = exp(x.v) * x.d;
		r.d = r.v * x.d;

		return r;
	}

	friend autodif log (const autodif& x) {
		autodif r;

		using std::log;
		r.v = log(x.v);
		r.d = x.d / x.v;

		return r;
	}

	friend autodif sqrt (const autodif& x) {
		autodif r;

		using std::sqrt;
		r.v = sqrt(x.v);
		// r.d = 1./(2. * sqrt(x.v)) * x.d;
		r.d = x.d / (2. * r.v);

		return r;
	}

	friend autodif sin (const autodif& x) {
		autodif r;

		using std::sin;
		using std::cos;
		r.v = sin(x.v);
		r.d = cos(x.v) * x.d;

		return r;
	}


	friend autodif cos (const autodif& x) {
		autodif r;

		using std::sin;
		using std::cos;
		r.v = cos(x.v);
		r.d = -sin(x.v) * x.d;

		return r;
	}

	friend autodif tan (const autodif& x) {
		autodif r;
		T tmp;

		using std::tan;
		using std::cos;
		r.v = tan(x.v);
		tmp = cos(x.v);
		tmp = 1. / (tmp * tmp);
		r.d = tmp * x.d;

		return r;
	}

	friend autodif asin (const autodif& x) {
		autodif r, tmp;

		using std::asin;
		using std::sqrt;
		r.v = asin(x.v);
		r.d = (1. / sqrt(1. - x.v * x.v)) * x.d;

		return r;
	}

	friend autodif acos (const autodif& x) {
		autodif r;

		using std::acos;
		using std::sqrt;
		r.v = acos(x.v);
		r.d = (-1. / sqrt(1. - x.v * x.v)) * x.d;

		return r;
	}

	friend autodif atan (const autodif& x) {
		autodif r;

		using std::atan;
		r.v = atan(x.v);
		r.d = (1. / (1. + x.v * x.v)) * x.d;

		return r;
	}

	friend autodif sinh (const autodif& x) {
		autodif r;

		using std::sinh;
		using std::cosh;
		r.v = sinh(x.v);
		r.d = cosh(x.v) * x.d;

		return r;
	}

	friend autodif cosh (const autodif& x) {
		autodif r;

		using std::sinh;
		using std::cosh;
		r.v = cosh(x.v);
		r.d = sinh(x.v) * x.d;

		return r;
	}

	friend autodif tanh (const autodif& x) {
		autodif r;
		T tmp;

		using std::tanh;
		using std::cosh;
		r.v = tanh(x.v);
		tmp = cosh(x.v);
		tmp = 1. / (tmp * tmp);
		r.d = tmp * x.d;

		return r;
	}

	friend autodif asinh (const autodif& x) {
		autodif r, tmp;

		// using std::asinh;
		using std::sqrt;
		r.v = asinh(x.v);
		r.d = (1. / sqrt(x.v * x.v + 1.)) * x.d;

		return r;
	}

	friend autodif acosh (const autodif& x) {
		autodif r;

		// using std::acosh;
		using std::sqrt;
		r.v = acosh(x.v);
		r.d = (1. / sqrt(x.v * x.v - 1.)) * x.d;

		return r;
	}

	friend autodif atanh (const autodif& x) {
		autodif r;

		// using std::atanh;
		r.v = atanh(x.v);
		r.d = (1. / (1. - x.v * x.v)) * x.d;

		return r;
	}

	// n-dimensional version
	static ub::vector<autodif> init (const ub::vector<T>& in) {
		int i, j;
		int n = in.size();
		ub::vector<autodif> out(n);

		for (i=0; i<n; i++) {
			out(i).v = in(i);
			out(i).d.resize(n);
			for (j=0; j<n; j++) {
				out(i).d(j) = (i==j) ? 1. : 0.;
			}
		}

		return out;
	}

	// 1-dimensional version
	static autodif init (const T& in) {
		int j;
		autodif out;

		out.v = in;
		out.d.resize(1);
		out.d(0) = 1.;

		return out;
	}

	// for functions R^n -> R^m
	static void split (const ub::vector<autodif>& in, ub::vector<T>& v, ub::matrix<T>& d) {
		int i, j, m, tmp;
		int n = in.size();

		m = in(0).d.size();
		for (i=1; i<n; i++) {
			tmp = in(i).d.size();
			if (tmp > m) m = tmp;
		}

		v.resize(n);
		d.resize(n, m);
		for (i=0; i<n; i++) {
			v(i) = in(i).v;
			tmp = in(i).d.size();
			for (j=0; j<tmp; j++) {
				d(i, j) = in(i).d(j);
			}
			for (j=tmp; j<m; j++) {
				d(i, j) = 0.;
			}
		}
	}

	// for functions R^n -> R
	static void split (const autodif& in, T& v, ub::vector<T>& d) {
		int j, m;

		m = in.d.size();
		d.resize(m);

		v = in.v;
		for (j=0; j<m; j++) {
			d(j) = in.d(j);
		}
	}

	static ub::vector<autodif>
	compress (const ub::vector<autodif>& in, ub::matrix<T>& save) {
		ub::vector<autodif> out;
		int i, j, m, tmp;
		int n = in.size();

		out.resize(n);

		m = in(0).d.size();
		for (i=1; i<n; i++) {
			tmp = in(i).d.size();
			if (tmp > m) m = tmp;
		}

		save.resize(n, m);

		for (i=0; i<n; i++) {
			out(i).v = in(i).v;
			tmp = in(i).d.size();
			for (j=0; j<tmp; j++) {
				save(i, j) = in(i).d(j);
			}
			for (j=tmp; j<m; j++) {
				save(i, j) = 0.;
			}
			out(i).d.resize(n);
			for (j=0; j<n; j++) {
				out(i).d(j) = (i==j) ? 1 : 0;
			}
		}

		return out;
	}

	static autodif
	expand (const autodif& in, const ub::matrix<T>& save) {
		autodif out;

		out.v = in.v;
		out.d = prod(in.d, save);

		return out;
	}

	static ub::vector<autodif>
	expand (const ub::vector<autodif>& in, const ub::matrix<T>& save) {
		ub::vector<autodif> out;
		int i;
		int s = in.size();

		out.resize(s);

		for (i=0; i<s; i++) {
			out(i) = expand(in(i), save);
		}

		return out;
	}
};

} // namespace kv

#endif //AUTODIF_HPP
