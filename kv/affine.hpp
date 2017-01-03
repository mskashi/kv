/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef AFFINE_HPP
#define AFFINE_HPP

// Affine Arithmetic

#include <iostream>
#include <stdexcept>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include <algorithm>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <kv/convert.hpp>


/*
 * define simplicity of affine arithmetic
 *
 *   value of AFFINE_SIMPLE                   | 0 | 1 | 2
 *  ------------------------------------------+---+---+---
 *   add dummy epsilin on linear operation    | o | x | x
 *   add dummy epsilin on nonlinear operation | o | o | x
 *
 *  default: 1
 */

#ifndef AFFINE_SIMPLE
#define AFFINE_SIMPLE 1
#endif

/*
 * select the method for multiplication
 *
 *  0: Stolfi's simple method (default, O(n))
 *  1: better multiplication (give smaller extra epsilon but slow, O(n^2))
 *  2: best multiplication (give smallest extra epsilon but slow, O(n^2))
 */

#ifndef AFFINE_MULT
#define AFFINE_MULT 0
#endif


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T> class affine;

template <class C, class T> struct acceptable_s<C, affine<T> > {
	static const bool value = boost::is_same<C, interval<T> >::value || boost::is_convertible<C, std::string>::value;
};

template <class C, class T> struct acceptable_n<C, affine<T> > {
	static const bool value = convertible<C, T>::value && (!acceptable_s<C, affine<T> >::value);
};

template <class C, class T> struct convertible<C, affine<T> > {
	static const bool value = acceptable_n<C, affine<T> >::value || acceptable_s<C, affine<T> >::value || boost::is_same<C, affine<T> >::value;
};

#if 0
template <class C, class T> struct acceptable_n<C, affine<T> > {
	static const bool value = convertible<C, T>::value && (!boost::is_same<C, interval<T> >::value) && (!boost::is_convertible<C, std::string>::value);
};

template <class C, class T> struct convertible<C, affine<T> > {
	static const bool value = convertible<C, T>::value || boost::is_same<C, interval<T> >::value || boost::is_convertible<C, std::string>::value || boost::is_same<C, affine<T> >::value;
};
#endif



template <class T> class affine {
	public:
	ub::vector<T> a;
	#if AFFINE_SIMPLE >= 1
	T er;
	#endif

	typedef T base_type;

	static int& maxnum() {
		static int m = 0;
		#pragma omp threadprivate(m)
		return m;
	}

	friend inline T rad(const affine& x) {
		int i, xs;
		T r(0.);
		T tmp;

		xs = x.a.size();

		rop<T>::begin();
		for (i=1; i<xs; i++) {
			// r = rop<T>::add_up(r, abs(x.a(i)));
			tmp = (x.a(i) >= 0.) ? x.a(i) : -x.a(i);
			r = rop<T>::add_up(r, tmp);
		}
		#if AFFINE_SIMPLE >= 1
		r = rop<T>::add_up(r, x.er);
		#endif
		rop<T>::end();

		return r;
	}

	friend inline interval<T> to_interval(const affine& x) {
		T t1, t2, t3;

		t1 = rad(x);
		rop<T>::begin();
		t2 = rop<T>::sub_down(x.a(0), t1);
		t3 = rop<T>::add_up(x.a(0), t1);
		rop<T>::end();
		return interval<T>(t2, t3);
	}


	affine() {
	}

	template <class C> explicit affine(const C& x, typename boost::enable_if_c< acceptable_n<C, affine>::value >::type* =0) {
		a.resize(1);
		a(0) = x;
		#if AFFINE_SIMPLE >= 1
		er = 0.;
		#endif
	}

	template <class C> explicit affine(const C& x, typename boost::enable_if_c< acceptable_s<C, affine>::value >::type* =0) {
		int i;
		interval<T> I(x);
		maxnum()++;
		a.resize(maxnum()+1);


		rop<T>::begin();
		a(0) = rop<T>::mul_up(rop<T>::add_up(I.upper(), I.lower()), T(0.5));
		a(maxnum()) = rop<T>::sub_up(a(0), I.lower());
		rop<T>::end();

		for (i=1; i<maxnum(); i++) a(i) = 0.;

		#if AFFINE_SIMPLE >= 1
		er = 0.;
		#endif
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, affine>::value, affine& >::type operator=(const C& x) {
		a.resize(1);
		a(0) = x;
		#if AFFINE_SIMPLE >= 1
		er = 0.;
		#endif

		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_s<C, affine>::value, affine& >::type operator=(const C& x) {
		int i;
		interval<T> I(x);
		maxnum()++;
		a.resize(maxnum()+1);

		rop<T>::begin();
		a(0) = rop<T>::mul_up(rop<T>::add_up(I.upper(), I.lower()), T(0.5));
		a(maxnum()) = rop<T>::sub_up(a(0), I.lower());
		rop<T>::end();

		for (i=1; i<maxnum(); i++) a(i) = 0.;

		#if AFFINE_SIMPLE >= 1
		er = 0.;
		#endif

		return *this;
	}

	T get_coef (int i) const {
		if (i >= a.size()) return T(0.);
		else return a(i);
	}

	T get_mid() const {
		return a(0);
	}

	T get_err() const {
		#if AFFINE_SIMPLE >= 1
		return er;
		#else
		return T(0.);
		#endif
	}

	friend affine operator+(const affine& x, const affine& y) {
		affine r;
		int xs, ys, i;
		T err(0.);

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		xs = x.a.size();
		ys = y.a.size();
		if (xs > ys) {
			#if AFFINE_SIMPLE >= 1
			r.a.resize(xs);
			#endif
			rop<T>::begin();
			for (i=0; i<ys; i++) {
				r.a(i) = rop<T>::add_down(x.a(i), y.a(i));
			}
			for (i=0; i<ys; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::add_up(x.a(i), y.a(i)), r.a(i)));
			}
			rop<T>::end();

			for (i=ys; i<xs; i++) {
				r.a(i) = x.a(i);
			}
			#if AFFINE_SIMPLE == 0
			for (i=xs; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		} else {
			#if AFFINE_SIMPLE >= 1
			r.a.resize(ys);
			#endif
			rop<T>::begin();
			for (i=0; i<xs; i++) {
				r.a(i) = rop<T>::add_down(x.a(i), y.a(i));
			}
			for (i=0; i<xs; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::add_up(x.a(i), y.a(i)), r.a(i)));
			}
			rop<T>::end();

			for (i=xs; i<ys; i++) {
				r.a(i) = y.a(i);
			}
			#if AFFINE_SIMPLE == 0
			for (i=ys; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		}
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(rop<T>::add_up(x.er, y.er), err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	// same as operator+, but do not add extra epsilon.
	// This function can be used for adding affine variables
	// which have no common epsilons.

	friend affine append(const affine& x, const affine& y) {
		affine r;
		int xs, ys, i;

		xs = x.a.size();
		ys = y.a.size();

		if (xs > ys) {
			r.a.resize(xs);
			for (i=0; i<ys; i++) {
				r.a(i) = x.a(i) + y.a(i);
			}
			for (i=ys; i<xs; i++) {
				r.a(i) = x.a(i);
			}
		} else {
			r.a.resize(ys);
			for (i=0; i<xs; i++) {
				r.a(i) = x.a(i) + y.a(i);
			}
			for (i=xs; i<ys; i++) {
				r.a(i) = y.a(i);
			}
		}
		#if AFFINE_SIMPLE >= 1
		r.er = x.er + y.er;
		#endif
		return r;
	}


	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator+(const affine& x, const C& y) {
		affine r;
		int xs, i;
		T err;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(xs);
		#endif

		rop<T>::begin();
		r.a(0) = rop<T>::add_down(x.a(0), (T)y);
		err = rop<T>::sub_up(rop<T>::add_up(x.a(0), (T)y), r.a(0));
		rop<T>::end();

		for (i=1; i<xs; i++) {
			r.a(i) = x.a(i);
		}
		#if AFFINE_SIMPLE == 0
		for (i=xs; i<maxnum(); i++) {
			r.a(i) = 0.;
		}
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(x.er, err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator+(const C& x, const affine& y) {
		affine r;
		int ys, i;
		T err;

		ys = y.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(ys);
		#endif

		rop<T>::begin();
		r.a(0) = rop<T>::add_down((T)x, y.a(0));
		err = rop<T>::sub_up(rop<T>::add_up((T)x, y.a(0)), r.a(0));
		rop<T>::end();

		for (i=1; i<ys; i++) {
			r.a(i) = y.a(i);
		}
		#if AFFINE_SIMPLE == 0
		for (i=ys; i<maxnum(); i++) {
			r.a(i) = 0.;
		}
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(y.er, err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator+(const affine& x, const C& y) {
		return x + affine(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator+(const C& x, const affine& y) {
		return affine(x) + y;
	}

	friend affine& operator+=(affine& x, const affine& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine& >::type operator+=(affine& x, const C& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine& >::type operator+=(affine& x, const C& y) {
		x = x + y;
		return x;
	}

	friend affine operator-(const affine& x, const affine& y) {
		affine r;
		int xs, ys, i;
		T err(0.);

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		xs = x.a.size();
		ys = y.a.size();
		if (xs > ys) {
			#if AFFINE_SIMPLE >= 1
			r.a.resize(xs);
			#endif
			rop<T>::begin();
			for (i=0; i<ys; i++) {
				r.a(i) = rop<T>::sub_down(x.a(i), y.a(i));
			}
			for (i=0; i<ys; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::sub_up(x.a(i), y.a(i)), r.a(i)));
			}
			rop<T>::end();

			for (i=ys; i<xs; i++) {
				r.a(i) = x.a(i);
			}
			#if AFFINE_SIMPLE == 0
			for (i=xs; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		} else {
			#if AFFINE_SIMPLE >= 1
			r.a.resize(ys);
			#endif
			rop<T>::begin();
			for (i=0; i<xs; i++) {
				r.a(i) = rop<T>::sub_down(x.a(i), y.a(i));
			}
			for (i=0; i<xs; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::sub_up(x.a(i), y.a(i)), r.a(i)));
			}
			rop<T>::end();

			for (i=xs; i<ys; i++) {
				r.a(i) = - y.a(i);
			}
			#if AFFINE_SIMPLE == 0
			for (i=ys; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		}
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(rop<T>::add_up(x.er, y.er), err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator-(const affine& x, const C& y) {
		affine r;
		int xs, i;
		T err;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(xs);
		#endif

		rop<T>::begin();
		r.a(0) = rop<T>::sub_down(x.a(0), (T)y);
		err = rop<T>::sub_up(rop<T>::sub_up(x.a(0), (T)y), r.a(0));
		rop<T>::end();

		for (i=1; i<xs; i++) {
			r.a(i) = x.a(i);
		}
		#if AFFINE_SIMPLE == 0
		for (i=xs; i<maxnum(); i++) {
			r.a(i) = 0.;
		}
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(x.er, err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator-(const C& x, const affine& y) {
		affine r;
		int ys, i;
		T err;

		ys = y.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(ys);
		#endif

		rop<T>::begin();
		r.a(0) = rop<T>::sub_down((T)x, y.a(0));
		err = rop<T>::sub_up(rop<T>::sub_up((T)x, y.a(0)), r.a(0));
		rop<T>::end();

		for (i=1; i<ys; i++) {
			r.a(i) = - y.a(i);
		}
		#if AFFINE_SIMPLE == 0
		for (i=ys; i<maxnum(); i++) {
			r.a(i) = 0.;
		}
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(y.er, err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator-(const affine& x, const C& y) {
		return x - affine(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator-(const C& x, const affine& y) {
		return affine(x) - y;
	}

	friend affine& operator-=(affine& x, const affine& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine& >::type operator-=(affine& x, const C& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine& >::type operator-=(affine& x, const C& y) {
		x = x - y;
		return x;
	}

	friend affine operator-(const affine& x) {
		affine r;

		r.a = - x.a;
		#if AFFINE_SIMPLE >= 1
		r.er = x.er;
		#endif
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator*(const affine& x, const C& y) {
		affine r;
		int xs, i;
		T err(0.);

		xs = x.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(xs);
		#endif

		rop<T>::begin();
		for (i=0; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), (T)y);
		}
		for (i=0; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), (T)y), r.a(i)));
		}
		rop<T>::end();

		#if AFFINE_SIMPLE == 0
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(rop<T>::mul_up(x.er, (T)((y >= 0.) ? y : -y)), err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator*(const C& x, const affine& y) {
		affine r;
		int ys, i;
		T err(0.);

		ys = y.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(ys);
		#endif

		rop<T>::begin();
		for (i=0; i<ys; i++) {
			r.a(i) = rop<T>::mul_down((T)x, y.a(i));
		}
		for (i=0; i<ys; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up((T)x, y.a(i)), r.a(i)));
		}
		rop<T>::end();

		#if AFFINE_SIMPLE == 0
		for (i=ys; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		r.er = rop<T>::add_up(rop<T>::mul_up(y.er, (T)((x >= 0.) ? x : -x)), err);
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator*(const affine& x, const C& y) {
		return x * affine(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator*(const C& x, const affine& y) {
		return affine(x) * y;
	}


	#if AFFINE_MULT >= 1

	static interval<T> bestmult_error(const affine& x, const affine& y)
	{
		int n, i, j, s;
		bool f1, f2;
		ub::vector<T> vx, vy;
		kv::interval<T> tmp, C, D, E, X, Y;

		n = std::max(x.a.size(), y.a.size()) - 1;

		#if AFFINE_SIMPLE >= 1
		vx.resize(n+3);
		vy.resize(n+3);
		#else
		vx.resize(n+1);
		vy.resize(n+1);
		#endif

		for (i=0; i<x.a.size(); i++) {
			vx(i) = x.a(i);
		}
		for (i=x.a.size(); i<=n; i++) {
			vx(i) = 0.;
		}
		for (i=0; i<y.a.size(); i++) {
			vy(i) = y.a(i);
		}
		for (i=y.a.size(); i<=n; i++) {
			vy(i) = 0.;
		}

		#if AFFINE_SIMPLE >= 1
		vx(n+1) = x.er;
		vy(n+1) = 0.;
		vx(n+2) = 0.;
		vy(n+2) = y.er;
		n += 2;
		#endif

		E = 0.;

		#if AFFINE_MULT == 1
		for (i=1; i<=n; i++) {
			X = vx(i);
			Y = vy(i);
			E += X * Y * interval<T>(0., 1.);
			for (j=i+1; j<=n; j++) {
				tmp = X * vy(j) + Y * vx(j);
				E += tmp * interval<T>(-1., 1.);
			}
		}
		#else // AFFINE_MULT == 1

		for (i=1; i<=n; i++) {
			X = vx(i);
			Y = vy(i);
			if (X == 0. && Y == 0.) continue;
			C = D = 0.;
			for (j=1; j<=n; j++) {
				if (j == i) continue;
				if (vy(j) == 0. && vx(j) == 0.) continue;
				tmp = X * vy(j) - Y * vx(j);
				if (tmp.lower() >= 0.) {
					C += vx(j);
					D += vy(j);
				} else if (tmp.upper() <= 0.) {
					C -= vx(j);
					D -= vy(j);
				} else {
					if ((X.upper() > 0. && vx(j) < 0.)
					|| (X.lower() < 0. && vx(j) > 0.)) {
						X -= vx(j);
						Y -= vy(j);
					} else {
						X += vx(j);
						Y += vy(j);
					}
				}
			}
			E = interval<T>::hull(E, X * Y + (X * D + Y * C) + C * D);
			E = interval<T>::hull(E, X * Y - (X * D + Y * C) + C * D);
			if (!zero_in(X * Y) ){
				tmp = -0.5 * (X * D + Y * C) / (X * Y);
				if (overlap(tmp, interval<T>(-1., 1.))) {
					E = interval<T>::hull(E, X * Y * tmp * tmp + (X * D + Y * C) * tmp + C * D);
				}
			}
		}

		#endif // AFFINE_MULT == 1

		return E;
	}

	#endif // AFFINE_MULT >= 1

	friend affine operator*(const affine& x, const affine& y) {
		affine r;
		int i, j, xs, ys;
		T err;
		T tmp_u, tmp_l;

		// if (&x == &y) return square(x);

		xs = x.a.size();
		ys = y.a.size();

		#if AFFINE_SIMPLE != 2
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(std::max(xs, ys));
		#endif

		rop<T>::begin();
		r.a(0) = rop<T>::mul_down(x.a(0), y.a(0));
		err = rop<T>::sub_up(rop<T>::mul_up(x.a(0), y.a(0)), r.a(0));
		rop<T>::end();

		if (xs > ys) {
			rop<T>::begin();
			for (i=1; i<ys; i++) {
				r.a(i) = rop<T>::add_down(rop<T>::mul_down(y.a(0), x.a(i)), rop<T>::mul_down(x.a(0), y.a(i)));
			}
			for (i=ys; i<xs; i++) {
				r.a(i) = rop<T>::mul_down(y.a(0), x.a(i));
			}
			for (i=1; i<ys; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::add_up(rop<T>::mul_up(y.a(0), x.a(i)), rop<T>::mul_up(x.a(0), y.a(i))), r.a(i)));
			}
			for (i=ys; i<xs; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(y.a(0), x.a(i)), r.a(i)));
			}
			rop<T>::end();
			#if AFFINE_SIMPLE != 2
			for (i=xs; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		} else {
			rop<T>::begin();
			for (i=1; i<xs; i++) {
				r.a(i) = rop<T>::add_down(rop<T>::mul_down(y.a(0), x.a(i)), rop<T>::mul_down(x.a(0), y.a(i)));
			}
			for (i=xs; i<ys; i++) {
				r.a(i) = rop<T>::mul_down(x.a(0), y.a(i));
			}
			for (i=1; i<xs; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::add_up(rop<T>::mul_up(y.a(0), x.a(i)), rop<T>::mul_up(x.a(0), y.a(i))), r.a(i)));
			}
			for (i=xs; i<ys; i++) {
				err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(0), y.a(i)), r.a(i)));
			}
			rop<T>::end();
			#if AFFINE_SIMPLE != 2
			for (i=ys; i<maxnum(); i++) {
				r.a(i) = 0.;
			}
			#endif
		}

		#if AFFINE_MULT >= 1

		interval<T> E;
		E = r.a(0) + bestmult_error(x, y);
		rop<T>::begin();
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(E.upper(), E.lower()), 0.5);
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), E.lower()));
		rop<T>::end();

		#else // AFFINE_MULT >= 1

		tmp_l = rad(x);
		tmp_u = rad(y);
		rop<T>::begin();
		err = rop<T>::add_up(err, rop<T>::mul_up(tmp_l, tmp_u));
		rop<T>::end();

		#endif // AFFINE_MULT >= 1

		#if AFFINE_SIMPLE >= 1
		rop<T>::begin();
		err = rop<T>::add_up(err, rop<T>::add_up(rop<T>::mul_up((y.a(0) >= 0.) ? y.a(0) : -y.a(0), x.er), rop<T>::mul_up((x.a(0) >= 0.) ? x.a(0) : -x.a(0), y.er)));
		rop<T>::end();
		#endif

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine& operator*=(affine& x, const affine& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine& >::type operator*=(affine& x, const C& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine& >::type operator*=(affine& x, const C& y) {
		x = x * y;
		return x;
	}


	friend affine inv(const affine& x) {
		affine r;
		T err;
		interval<T> I, tmp, range;
		T a, b, l, u;
		int i, xs;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		I = to_interval(x);
		l = I.lower();
		u = I.upper();
		if (l < 0. && u > 0.) {
			throw std::domain_error("affine: division by 0");
		}

		a = -1. /(l * u);
		tmp = a; // tmp is used to force interval calculation
		if (u > 0.) {
			range = 2. * sqrt(-tmp);
		} else {
			range = -2. * sqrt(-tmp);
		}
		tmp = l;
		range = interval<T>::hull(range, 1./tmp - a * tmp);
		tmp = u;
		range = interval<T>::hull(range, 1./tmp - a * tmp);

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine operator/(const affine& x, const affine& y) {
		return x * inv(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator/(const C& x, const affine& y) {
		return x * inv(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type operator/(const affine& x, const C& y) {
		affine r;
		int xs, i;
		T err(0.);

		xs = x.a.size();

		#if AFFINE_SIMPLE == 0
		maxnum()++;
		r.a.resize(maxnum()+1);
		#else
		r.a.resize(xs);
		#endif

		rop<T>::begin();
		for (i=0; i<xs; i++) {
			r.a(i) = rop<T>::div_down(x.a(i), (T)y);
		}
		for (i=0; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::div_up(x.a(i), (T)y), r.a(i)));
		}
		rop<T>::end();

		#if AFFINE_SIMPLE == 0
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		#if AFFINE_SIMPLE >= 1
		// r.er = x.er / abs(y) + err;
		rop<T>::begin();
		r.er = rop<T>::add_up(err, rop<T>::div_up(x.er, (T)((y >= 0.) ? y : -y)));
		rop<T>::end();
		#else
		r.a(maxnum()) = err;
		#endif

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator/(const affine& x, const C& y) {
		return x / affine(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type operator/(const C& x, const affine& y) {
		return affine(x) / y;
	}

	friend affine& operator/=(affine& x, const affine& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine& >::type operator/=(affine& x, const C& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine& >::type operator/=(affine& x, const C& y) {
		x = x / y;
		return x;
	}


	friend affine sqrt(const affine& x) {
		affine r;
		T err;
		interval<T> I, tmp, range;
		T a, b, l, u;
		int i, xs;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		I = to_interval(x);
		l = I.lower();
		u = I.upper();
		if (l < 0.) {
			throw std::domain_error("affine: sqrt of negative value");

		}

		using std::sqrt;
		a = 1. /(sqrt(l) + sqrt(u));
		tmp = a;
		range = 1. / (4. * tmp);
		tmp = l;
		range = interval<T>::hull(range, sqrt(tmp) - a * tmp);
		tmp = u;
		range = interval<T>::hull(range, sqrt(tmp) - a * tmp);

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine square(const affine& x) {
		affine r;
		T err;
		interval<T> I, tmp, range;
		T a, b, l, u;
		int i, xs;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		I = to_interval(x);
		l = I.lower();
		u = I.upper();

		a = l + u;
		tmp = a;
		range = - tmp * tmp * 0.25;
		tmp = l;
		// range = hull(range, tmp * tmp - a * tmp);
		range = interval<T>::hull(range, tmp * (tmp - a));
		tmp = u;
		range = interval<T>::hull(range, tmp * (tmp - a));

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine exp(const affine& x) {
		affine r;
		T err;
		interval<T> I, tmp, range;
		T a, b, l, u;
		int i, xs;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		I = to_interval(x);
		l = I.lower();
		u = I.upper();

		using std::exp;
		if (u == l) a = exp(u);
		else a = (exp(u) - exp(l)) / (u - l);
		tmp = a;
		range = tmp * (1. - log(tmp));
		tmp = l;
		range = interval<T>::hull(range, exp(tmp) - a * tmp);
		tmp = u;
		range = interval<T>::hull(range, exp(tmp) - a * tmp);

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine log(const affine& x) {
		affine r;
		T err;
		interval<T> I, tmp, range;
		T a, b, l, u;
		int i, xs;

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		I = to_interval(x);
		l = I.lower();
		u = I.upper();

		if (l <= 0.) {
			throw std::domain_error("affine: log of nagative value");
		}

		using std::log;
		if (u == l) a = 1. / u;
		else a = (log(u) - log(l)) / (u - l);
		tmp = a;
		range = log(1. / tmp) - 1.;
		tmp = l;
		range = interval<T>::hull(range, log(tmp) - a * tmp);
		tmp = u;
		range = interval<T>::hull(range, log(tmp) - a * tmp);

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	friend affine abs(const affine& x) {
		affine r;
		T err;
		interval<T> I, range;
		T a, b, l, u;
		int i, xs;

		I = to_interval(x);
		l = I.lower();
		u = I.upper();

		if (l >= 0.) {
			return x;
		}
		if (u <= 0.) {
			return -x;
		}

		xs = x.a.size();

		#if AFFINE_SIMPLE == 2
		r.a.resize(xs);
		#else
		maxnum()++;
		r.a.resize(maxnum()+1);
		#endif

		a = (u + l) / (u - l);

		range = 0.;
		range = interval<T>::hull(range, u - a * u);
		range = interval<T>::hull(range, -l - a * l);

		rop<T>::begin();
		b = rop<T>::mul_up(rop<T>::add_up(range.upper(), range.lower()), T(0.5));
		err = rop<T>::sub_up(b, range.lower());
		for (i=1; i<xs; i++) {
			r.a(i) = rop<T>::mul_down(x.a(i), a);
		}
		for (i=1; i<xs; i++) {
			err = rop<T>::add_up(err, rop<T>::sub_up(rop<T>::mul_up(x.a(i), a), r.a(i)));
		}
		#if AFFINE_SIMPLE != 2
		for (i=xs; i<maxnum(); i++) r.a(i) = 0.;
		#endif
		l = rop<T>::add_down(rop<T>::mul_down(x.a(0), a), b);
		u = rop<T>::add_up(rop<T>::mul_up(x.a(0), a), b);
		r.a(0) = rop<T>::mul_up(rop<T>::add_up(l, u), T(0.5));
		err = rop<T>::add_up(err, rop<T>::sub_up(r.a(0), l));
		#if AFFINE_SIMPLE >= 1
		// err += abs(a) * x.er;
		err = rop<T>::add_up(err, rop<T>::mul_up(((a >= 0.) ? a : -a), x.er));
		#endif
		rop<T>::end();

		#if AFFINE_SIMPLE == 2
		r.er = err;
		#else
		r.a(maxnum()) = err;
		# if AFFINE_SIMPLE == 1
		r.er = 0.;
		# endif
		#endif

		return r;
	}

	// lazy implementation of integer pow
	friend affine pow(const affine& x, int y) {
		affine r, xp;
		int a, tmp;

		if (y == 0) return affine(1.);

		a = (y >= 0) ? y : -y;

		tmp = a;
		r = 1.;
		xp = x;
		while (tmp != 0) {
			if (tmp % 2 != 0) {
				r *= xp;
			}
			tmp /= 2;
			xp = xp * xp;
		}

		if (y < 0) {
			r = 1. / r;
		}

		return r;
	}

	friend affine pow(const affine& x, const affine& y) {
		return exp(y * log(x));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value && ! boost::is_integral<C>::value, affine >::type pow(const affine& x, const C& y) {
		return pow(x, affine(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type pow(const affine& x, const C& y) {
		return pow(x, affine(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, affine>::value, affine >::type pow(const C& x, const affine& y) {
		return pow(affine(x), y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, affine>::value, affine >::type pow(const C& x, const affine& y) {
		return pow(affine(x), y);
	}

	// lazy implementation of sin
	friend affine sin(const affine& x) {
		// return sin((interval<T>)(x.a(0))) + cos(to_interval(x)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = sin((interval<T>)(x.a(0)));
		r2 = cos(to_interval(x));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of cos
	friend affine cos(const affine& x) {
		// return cos((interval<T>)(x.a(0))) - sin(to_interval(x)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = cos((interval<T>)(x.a(0)));
		r2 = -sin(to_interval(x));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of tan
	friend affine tan(const affine& x) {
		// return tan((interval<T>)(x.a(0))) + pow(cos(to_interval(x)), -2) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = tan((interval<T>)(x.a(0)));
		r2 = pow(cos(to_interval(x)), -2);
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of asin
	friend affine asin(const affine& x) {
		// return asin((interval<T>)(x.a(0))) + (1. / sqrt(1. - pow(to_interval(x), 2))) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = asin((interval<T>)(x.a(0)));
		r2 = 1. / sqrt(1. - pow(to_interval(x), 2));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of acos
	friend affine acos(const affine& x) {
		// return acos((interval<T>)(x.a(0))) + (- 1. / sqrt(1. - pow(to_interval(x), 2))) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = acos((interval<T>)(x.a(0)));
		r2 = - 1. / sqrt(1. - pow(to_interval(x), 2));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of atan
	friend affine atan(const affine& x) {
		// return atan((interval<T>)(x.a(0))) + (1. / sqrt(1. + pow(to_interval(x), 2))) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = atan((interval<T>)(x.a(0)));
		r2 = 1. / sqrt(1. + pow(to_interval(x), 2));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of sinh
	friend affine sinh(const affine& x) {
		// return sinh((interval<T>)(x.a(0))) + cosh(to_interval(x)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = sinh((interval<T>)(x.a(0)));
		r2 = cosh(to_interval(x));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of cosh
	friend affine cosh(const affine& x) {
		// return cosh((interval<T>)(x.a(0))) + sinh(to_interval(x)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = cosh((interval<T>)(x.a(0)));
		r2 = sinh(to_interval(x));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of tanh
	friend affine tanh(const affine& x) {
		// return tanh((interval<T>)(x.a(0))) + pow(cosh(to_interval(x)), -2) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = tanh((interval<T>)(x.a(0)));
		r2 = pow(cosh(to_interval(x)), -2);
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of asinh
	friend affine asinh(const affine& x) {
		// return asinh((interval<T>)(x.a(0))) + (1 / sqrt(pow(to_interval(x), 2) + 1)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = asinh((interval<T>)(x.a(0)));
		r2 = 1 / sqrt(pow(to_interval(x), 2) + 1);
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of acosh
	friend affine acosh(const affine& x) {
		// return acosh((interval<T>)(x.a(0))) + (1 / sqrt(pow(to_interval(x), 2) - 1)) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = acosh((interval<T>)(x.a(0)));
		r2 = 1 / sqrt(pow(to_interval(x), 2) - 1);
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	// lazy implementation of atanh
	friend affine atanh(const affine& x) {
		// return atanh((interval<T>)(x.a(0))) + (1 / (1 - pow(to_interval(x), 2))) * (x - x.a(0));

		affine tmp;
		interval<T> r, r2;
		T m;

		tmp = x;
		tmp.a(0) = 0.;
		r = atanh((interval<T>)(x.a(0)));
		r2 = 1 / (1 - pow(to_interval(x), 2));
		m = mid(r2);
		r += (r2 - m) * to_interval(tmp);
		
		return r + m * tmp;
	}

	friend std::ostream& operator<<(std::ostream& s, const affine& x) {
		int i;

		s << "[(" << x.a(0) << ")";
		for (i=1; i<x.a.size(); i++) {
			s << "+(" << x.a(i) << ")e" << i;
		}
		#if AFFINE_SIMPLE >= 1
		s << "+(" << x.er << ")er";
		#endif
		s << "]";

		return s;
	}

	friend inline void split(const affine& x, int n, affine& y, affine& z) {
		int i;
		int s = x.a.size();
	
		y.a.resize(s);
		z.a.resize(s);
		for (i=0; i<=n; i++) {
			y.a(i) = x.a(i);
			z.a(i) = 0.;
		}
		for (i=n+1; i<s; i++) {
			y.a(i) = 0.;
			z.a(i) = x.a(i);
		}
		#if AFFINE_SIMPLE >= 1
		y.er = 0.;
		z.er = x.er;
		#endif
	}

	void resize() {
		int i;
		ub::vector<T> r;

		r.resize(maxnum()+1);
		for (i=0; i<=maxnum(); i++) {
			r(i) = a(i);
		}

		a = r;
	}
};


template <class T> inline ub::vector< interval<T> > to_interval(const ub::vector< affine<T> >& x) {
	int s = x.size();
	ub::vector< interval<T> > r;
	int i;

	r.resize(s);
	for (i=0; i<s; i++) r(i) = to_interval(x(i));

	return r;
}


/*
 * epsilon_reduce(x, n, n_limit)
 *  x: vector of affine (overwrited)
 *  reduce the number of epsilons to n if the number of epsilons > n_limit.
 *  
 *  s: size of x
 *  keep n-s "important" epsilons and convert other "non-important" epsilons
 *  to s-dimensional rectangular.
 */


// a class to store column vector and its "importance"

template <class T> class ep_reduce_v {
	public:
	ub::vector<T> v;
	T score;
	void calc_score() {
		int s = v.size();
		int i, j;
		T m1, m2, tmp;
		using std::abs;
		m1 = abs(v(0));
		m2 = abs(v(1));
		if (m2 > m1) {
			tmp = m2; m2 = m1; m1 = tmp;
		}
		for (i=2; i<s; i++) {
			tmp = abs(v(i));
			if (tmp > m1) {
				m2 = m1; m1 = tmp;
			} else if (tmp > m2) {
				m2 = tmp;
			}
		}
		if (m1 == 0.) score = 0.;
		else score = (m1*m2)/(m1+m2);
	}
};


// function to sort column vector by score

template <class T> inline bool ep_reduce_cmp(ep_reduce_v<T>* a, ep_reduce_v<T>* b) {
#ifdef EP_REDUCE_REVERSE
	return a->score < b->score;
#else
	return a->score > b->score;
#endif
}

template <class T> inline void epsilon_reduce(ub::vector< affine<T> >& x, int n, int n_limit = 0) {
	int s = x.size();
	int m = affine<T>::maxnum();
	int i, j;
	std::vector< ep_reduce_v<T> > a;
	std::vector< ep_reduce_v<T>* > pa;
	ub::vector< affine<T> > r;
	T tmp;

	if (n_limit < n) n_limit = n;

	if (m <= n_limit) return;
	if (n < s) return; // impossible

	a.resize(m);
	pa.resize(m);

	for (i=1; i<=m; i++) {
		a[i-1].v.resize(s);
		for (j=0; j<s; j++) {
			a[i-1].v(j) = (i < x(j).a.size()) ? x(j).a(i) : (T)0.;
		}
		a[i-1].calc_score();
		pa[i-1] = &(a[i-1]);
	}

#ifdef EP_REDUCE_REVERSE
	std::partial_sort(pa.begin(), pa.begin()+m-n+s, pa.end(), ep_reduce_cmp<T>);
#else
	std::partial_sort(pa.begin(), pa.begin()+n-s, pa.end(), ep_reduce_cmp<T>);
#endif

	r.resize(s);
	for (i=0; i<s; i++) {
		r(i).a.resize(n+1);
		r(i).a(0) = x(i).a(0);
		for (j=0; j<n-s; j++) {
#ifdef EP_REDUCE_REVERSE
			r(i).a(j+1) = pa[m-1-j]->v(i);
#else
			r(i).a(j+1) = pa[j]->v(i);
#endif
		}
		tmp = 0.;
		rop<T>::begin();
		for (j=n-s; j<m; j++) {
			using std::abs;
#ifdef EP_REDUCE_REVERSE
			tmp = rop<T>::add_up(tmp, abs(pa[m-1-j]->v(i)));
#else
			tmp = rop<T>::add_up(tmp, abs(pa[j]->v(i)));
#endif
		}
		#if AFFINE_SIMPLE >= 1
		tmp = rop<T>::add_up(tmp, x(i).er);
		#endif
		rop<T>::end();
		for (j=n-s; j<n; j++) {
			r(i).a(j+1) = 0.;
		}
		r(i).a(n-s+i+1) = tmp;
		#if AFFINE_SIMPLE >= 1
		r(i).er = 0.;
		#endif
	}

	x = r;
	affine<T>::maxnum() = n;
}


template <class T> struct constants< affine<T> > {
	static affine<T> pi() {
		static const affine<T> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static affine<T> e() {
		static const affine<T> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static affine<T> ln2() {
		static const affine<T> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}
	static affine<T> str(const std::string& s) {
		return affine<T>(s, s);
	}
	static affine<T> str(const std::string& s1, const std::string& s2) {
		return affine<T>(s1, s2);
	}
};


} // namespace kv

#endif //AFFINE_HPP
