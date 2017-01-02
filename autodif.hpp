#ifndef AUTODIF_HPP
#define AUTODIF_HPP

// Automatic Differentiation

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace ub = boost::numeric::ublas;


template <class T> class autodif {
	public:
	T v;
	ub::vector<T> d;

	autodif() {
	}

	// CからTへ暗黙の型変換が可能であるようなCに対してのみ生成される。
	template <class C> autodif(const C& x, typename boost::enable_if< boost::is_convertible<C, T> >::type* =0) {
		v = x;
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

	// CからTへ暗黙の型変換が可能であるようなCに対してのみ生成される。
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator+(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v + b;
		r.d = a.d;

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator+(const C& a, const autodif& b) {
		autodif r;

		r.v = a + b.v;
		r.d = b.d;

		return r;
	}

#if 0
	friend autodif& operator+=(autodif& a, const autodif& b) {
		a.v += b.v;

		if (a.d.size() == 0) {
			a.d = b.d;
		} else if (b.d.size() == 0) {
		} else {
			a.d += b.d;
		}

		return a;
	}
#endif

	// Cからautodifへ暗黙の型変換が可能であるようなCに対してのみ生成される。
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, autodif>, autodif& >::type operator+=(autodif& a, const C& b) {
		a = a + b;
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

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator-(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v - b;
		r.d = a.d;

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator-(const C& a, const autodif& b) {
		autodif r;

		r.v = a - b.v;
		r.d = - b.d;

		return r;
	}

#if 0
	friend autodif& operator-=(autodif& a, const autodif& b) {
		a.v -= b.v;

		if (a.d.size() == 0) {
			a.d = -b.d;
		} else if (b.d.size() == 0) {
		} else {
			a.d -= b.d;
		}

		return a;
	}
#endif

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, autodif>, autodif& >::type operator-=(autodif& a, const C& b) {
		a = a - b;
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

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator*(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v * b;
		r.d = b * a.d;

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator*(const C& a, const autodif& b) {
		autodif r;

		r.v = a * b.v;
		r.d = a * b.d;

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, autodif>, autodif& >::type operator*=(autodif& a, const C& b) {
		a = a * b;
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

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator/(const autodif& a, const C& b) {
		autodif r;

		r.v = a.v / b;
		r.d = a.d / b;

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, autodif >::type operator/(const C& a, const autodif& b) {
		autodif r;

		r.v = a / b.v;
		r.d = b.d * (-a/(b.v*b.v));

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, autodif>, autodif& >::type operator/=(autodif& a, const C& b) {
		a = a / b;
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

	friend autodif sin (const autodif& x) {
		autodif r;

		r.v = sin(x.v);
		r.d = cos(x.v) * x.d;

		return r;
	}

	friend autodif cos (const autodif& x) {
		autodif r;

		r.v = cos(x.v);
		r.d = -sin(x.v) * x.d;

		return r;
	}

	friend autodif exp (const autodif& x) {
		autodif r;

		r.v = exp(x.v);
		// r.d = exp(x.v) * x.d;
		r.d = r.v * x.d;

		return r;
	}

	friend autodif log (const autodif& x) {
		autodif r;

		r.v = log(x.v);
		r.d = x.d / x.v;

		return r;
	}

	friend autodif sqrt (const autodif& x) {
		autodif r;

		r.v = sqrt(x.v);
		// r.d = 1./(2. * sqrt(x.v)) * x.d;
		r.d = x.d / (2. * r.v);

		return r;
	}
};

// n-dimensional version
template<class T> inline ub::vector< autodif<T> > init_dif (const ub::vector<T>& in) {
 	int i, j;
	int n = in.size();
	ub::vector< autodif<T> > out(n);

	for (i=0; i<n; i++) {
		out(i).v = in(i);
		out(i).d.resize(n);
		for (j=0; j<n; j++) {
			out(i).d(j) = (i==j) ? 1 : 0;
		}
	}

	return out;
}

// 1-dimensional version
template<class T> inline autodif<T> init_dif (const T& in) {
 	int j;
	autodif<T> out;

	out.v = in;
	out.d.resize(1);
	out.d(0) = 1;

	return out;
}

template<class T> inline void split_dif (const ub::vector< autodif<T> >& in, ub::vector<T>& v, ub::matrix<T>& d) {
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

// R^n -> R なversion
template<class T> inline void split_dif (const autodif<T>& in, T& v, ub::vector<T>& d) {
 	int j, m;

	m = in.d.size();
	d.resize(m);

	v = in.v;
	for (j=0; j<m; j++) {
		d(j) = in.d(j);
	}
}

template<class T> inline void split_dif2 (const ub::vector< autodif<T> >& in, ub::vector<T>& v, ub::matrix<T>& d, int m) {
 	int i, j, tmp;
	int n = in.size();

	v.resize(n);
	d.resize(n, m);
	for (i=0; i<n; i++) {
		v(i) = in(i).v;
		tmp = in(i).d.size();
		if (m < tmp) tmp = m;
		for (j=0; j<tmp; j++) {
			d(i, j) = in(i).d(j);
		}
		for (j=tmp; j<m; j++) {
			d(i, j) = 0;
		}
	}
}

template<class T> inline ub::vector< autodif<T> >
compress_dif
(const ub::vector< autodif<T> >& in, ub::matrix<T>& save)
{
	ub::vector< autodif<T> > out;
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

template<class T> inline autodif<T>
expand_dif
(const autodif<T>& in, const ub::matrix<T>& save)
{
	autodif<T> out;

	out.v = in.v;
	out.d = prod(in.d, save);

	return out;
}

template<class T> inline ub::vector< autodif<T> >
expand_dif
(const ub::vector< autodif<T> >& in, const ub::matrix<T>& save)
{
	ub::vector< autodif<T> > out;
	int i;
	int s = in.size();

	out.resize(s);

	for (i=0; i<s; i++) {
		out(i) = expand_dif(in(i), save);
	}

	return out;
}

#if 0
template<class T> inline autodif<T> sin (const autodif<T>& x) {
	autodif<T> r;

	r.v = sin(x.v);
	r.d = cos(x.v) * x.d;

	return r;
}

template<class T> inline autodif<T> cos (const autodif<T>& x) {
	autodif<T> r;

	r.v = cos(x.v);
	r.d = -sin(x.v) * x.d;

	return r;
}
#endif

#endif //AUTODIF_HPP
