/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef PSA_HPP
#define PSA_HPP

// Power Series Arithmetic Type I and II with recording

#include <iostream>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/convert.hpp>

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T> class psa;

template <class C, class T> struct convertible<C, psa<T> > {
	static const bool value = convertible<C, T>::value || boost::is_same<C, psa<T> >::value;
};

template <class C, class T> struct acceptable_n<C, psa<T> > {
	static const bool value = convertible<C, T>::value;
};



template <class T> class psa {
	public:
	ub::vector<T> v;

	typedef T base_type;

	static int& mode() {
		static int m = 1;
		#pragma omp threadprivate (m)
		return m;
	}

	static T& domain() {
#ifdef _OPENMP // hack for non-POD thread local storage
		static T* d = NULL;
		#pragma omp threadprivate (d)
		if (d == NULL) {
			d = new T();
		}
		return *d;
#else
		static T d;
		return d;
#endif
	}


	static std::list<psa>& history() {
#ifdef _OPENMP // hack for non-POD thread local storage
		static std::list<psa>* hist = NULL;
		#pragma omp threadprivate (hist)
		if (hist == NULL) {
			hist = new std::list<psa>();
		}
		return *hist;
#else
		static std::list<psa> hist;
		return hist;
#endif
	}

	static bool& record_history() {
		static bool rh = false;
		#pragma omp threadprivate (rh)
		return rh;
	}

	static bool& use_history() {
		static bool uh = false;
		#pragma omp threadprivate (uh)
		return uh;
	}

	psa() {
		v.resize(1);
		v(0) = 0.;
	}

	template <class C> explicit psa(const C& x, typename boost::enable_if_c< acceptable_n<C, psa>::value >::type* =0) {
		v.resize(1);
		v(0) = x;
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, psa>::value, psa& >::type operator=(const C& x) {
		v.resize(1);
		v(0) = x;
		return *this;
	}

	friend psa operator+(const psa& a, const psa& b) {
		psa r;

		if (a.v.size() == 1) {
			r.v = b.v;
			r.v(0) += a.v(0);
		} else if (b.v.size() == 1) {
			r.v = a.v;
			r.v(0) += b.v(0);
		} else {
			if (use_history() == true) {
				r = history().front();
				int old_size = r.v.size();
				if (mode() == 2) {
					old_size = std::min(old_size, (int)a.v.size() - 1);
				}
				int i;
				r.v.resize(a.v.size(), true);
				for (i=old_size; i<a.v.size(); i++) {
					r.v(i) = a.v(i) + b.v(i);
				}
			} else {
				r.v = a.v + b.v;
			}
		}

		if (use_history() == true) {
			history().pop_front();
		}

		if (record_history() == true) {
			if (mode() == 1) {
				history().push_back(r);
			} else {
				psa r2 = r;
				r2.v.resize(r.v.size() - 1, true);
				history().push_back(r2);
			}
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator+(const psa& a, const C& b) {
		psa r;

		r.v = a.v;
		r.v(0) += b;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator+(const C& a, const psa& b) {
		psa r;

		r.v = b.v;
		r.v(0) += a;

		return r;
	}

	friend psa& operator+=(psa& a, const psa& b) {
		a = a + b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa& >::type operator+=(psa& a, const C& b) {
		a.v(0) += b;
		return a;
	}

	friend psa operator-(const psa& a, const psa& b) {
		psa r;

		if (a.v.size() == 1) {
			r.v = - b.v;
			r.v(0) += a.v(0);
		} else if (b.v.size() == 1) {
			r.v = a.v;
			r.v(0) -= b.v(0);
		} else {
			if (use_history() == true) {
				r = history().front();
				int old_size = r.v.size();
				if (mode() == 2) {
					old_size = std::min(old_size, (int)a.v.size() - 1);
				}
				int i;
				r.v.resize(a.v.size(), true);
				for (i=old_size; i<a.v.size(); i++) {
					r.v(i) = a.v(i) - b.v(i);
				}
			} else {
				r.v = a.v - b.v;
			}
		}

		if (use_history() == true) {
			history().pop_front();
		}

		if (record_history() == true) {
			if (mode() == 1) {
				history().push_back(r);
			} else {
				psa r2 = r;
				r2.v.resize(r.v.size() - 1, true);
				history().push_back(r2);
			}
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator-(const psa& a, const C& b) {
		psa r;

		r.v = a.v;
		r.v(0) -= b;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator-(const C& a, const psa& b) {
		psa r;

		r.v = - b.v;
		r.v(0) += a;

		return r;
	}

	friend psa& operator-=(psa& a, const psa& b) {
		a = a - b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa& >::type operator-=(psa& a, const C& b) {
		a.v(0) -= b;
		return a;
	}

	friend psa operator-(const psa& a) {
		psa r;

		r.v = - a.v;

		return r;
	}

	friend psa operator*(const psa& a, const psa& b) {
		psa r;
		int i, j, s;
		T sum;
		int old_size;

		if (a.v.size() == 1) {
			if (use_history() == true) {
				r = history().front();
				old_size = r.v.size();
				if (mode() == 2) {
					old_size = std::min(old_size, (int)b.v.size() - 1);
				}
				r.v.resize(b.v.size(), true);
				for (i=old_size; i<b.v.size(); i++) {
					r.v(i) = a.v(0) * b.v(i);
				}
			} else {
				r.v = a.v(0) * b.v;
			}
		} else if (b.v.size() == 1) {
			if (use_history() == true) {
				r = history().front();
				old_size = r.v.size();
				if (mode() == 2) {
					old_size = std::min(old_size, (int)a.v.size() - 1);
				}
				r.v.resize(a.v.size(), true);
				for (i=old_size; i<a.v.size(); i++) {
					r.v(i) = a.v(i) * b.v(0);
				}
			} else {
				r.v = a.v * b.v(0);
			}
		} else {
			s = a.v.size();
			if (use_history() == true) {
				r = history().front();
				old_size = r.v.size();
				r.v.resize(s, true);
			} else {
				old_size = 0;
				r.v.resize(s);
			}
			for (i=old_size; i<s; i++) {
				sum = 0.;
				for (j=0; j<=i; j++) {
					sum += a.v(j) * b.v(i-j);
				}
				r.v(i) = sum;
			}

			if (mode() == 2) {
				// tmpの中身もhistoryが利用できそうだが、保留
				ub::vector<T> tmp(s);
				tmp(0) = r.v(s-1);
				for (i=1; i<s; i++) {
					sum = 0.;
					for (j=i; j<s; j++) {
						sum += a.v(j) * b.v(i-j+s-1);
					}
					tmp(i) = sum;
				}
				r.v(s-1) = polyrange(tmp, 0, s-1, domain());
			}
		}

		if (use_history() == true) {
			history().pop_front();
		}

		if (record_history() == true) {
			if (mode() == 1) {
				history().push_back(r);
			} else {
				psa r2 = r;
				r2.v.resize(r.v.size() - 1, true);
				history().push_back(r2);
			}
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator*(const psa& a, const C& b) {
		psa r;

		if (use_history() == true) {
			r = history().front();
			int old_size = r.v.size();
			if (mode() == 2) {
				old_size = std::min(old_size, (int)a.v.size() - 1);
			}
			int i;
			r.v.resize(a.v.size(), true);
			for (i=old_size; i<a.v.size(); i++) {
				r.v(i) = a.v(i) * b;
			}
		} else {
			// r.v = a.v * b;
			r.v = a.v * T(b); // assist for VC++
		}

		if (use_history() == true) {
			history().pop_front();
		}

		if (record_history() == true) {
			if (mode() == 1) {
				history().push_back(r);
			} else {
				psa r2 = r;
				r2.v.resize(r.v.size() - 1, true);
				history().push_back(r2);
			}
		}

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator*(const C& a, const psa& b) {
		psa r;

		if (use_history() == true) {
			r = history().front();
			int old_size = r.v.size();
			if (mode() == 2) {
				old_size = std::min(old_size, (int)b.v.size() - 1);
			}
			int i;
			r.v.resize(b.v.size(), true);
			for (i=old_size; i<b.v.size(); i++) {
				r.v(i) = a * b.v(i);
			}
		} else {
			// r.v = a * b.v;
			r.v = T(a) * b.v; // assist for VC++
		}

		if (use_history() == true) {
			history().pop_front();
		}

		if (record_history() == true) {
			if (mode() == 1) {
				history().push_back(r);
			} else {
				psa r2 = r;
				r2.v.resize(r.v.size() - 1, true);
				history().push_back(r2);
			}
		}

		return r;
	}

	friend psa& operator*=(psa& a, const psa& b) {
		a = a * b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa& >::type operator*=(psa& a, const C& b) {
		a = a * b;
		return a;
	}

	friend psa operator/(const psa& a, const psa& b) {
		return a * inv(b);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator/(const psa& a, const C& b) {
		psa r;

		// r.v = a.v / b;
		r.v = a.v / T(b); // assist for VC++

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa >::type operator/(const C& a, const psa& b) {
		return a * inv(b);
	}

	friend psa& operator/=(psa& a, const psa& b) {
		a = a / b;
		return a;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, psa>::value, psa& >::type operator/=(psa& a, const C& b) {
		a = a / b;
		return a;
	}

	friend psa inv(const psa& x) {
		T a, xn, xn2, range;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		r = 1./a;
		hn = 1.;
		xn = 1./a;
		range = evalrange(x);
		xn2 = 1./range;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				xn2 = -xn2 / range;
				r += xn2 * hn;
			} else {
				hn *= h;
				xn = -xn / a;
				xn2 = -xn2 / range;
				r += xn * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa sin (const psa& x) {
		T a, xn, xn2, range, fact_n, table[4], tmp;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		table[0] = sin(a);
		table[1] = cos(a);
		table[2] = -table[0];
		table[3] = -table[1];

		r = table[0];
		hn = 1.;
		fact_n = 1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				fact_n /= (double)i;
				range = evalrange(x);
				switch (i%4) {
				case 0: tmp = sin(range); break;
				case 1: tmp = cos(range); break;
				case 2: tmp = -sin(range); break;
				case 3: tmp = -cos(range); break;
				}
				r += fact_n * tmp * hn;
			} else {
				hn *= h;
				fact_n /= (double)i;
				r += fact_n * table[i % 4] * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa cos (const psa& x) {
		T a, xn, xn2, range, fact_n, table[4], tmp;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		table[0] = cos(a);
		table[1] = -sin(a);
		table[2] = -table[0];
		table[3] = -table[1];

		r = table[0];
		hn = 1.;
		fact_n = 1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				fact_n /= (double)i;
				range = evalrange(x);
				switch (i%4) {
				case 0: tmp = cos(range); break;
				case 1: tmp = -sin(range); break;
				case 2: tmp = -cos(range); break;
				case 3: tmp = sin(range); break;
				}
				r += fact_n * tmp * hn;
			} else {
				hn *= h;
				fact_n /= (double)i;
				r += fact_n * table[i % 4] * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa exp (const psa& x) {
		T a, xn, xn2, range, fact_n, ea;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		ea = exp(a);

		r = ea;
		hn = 1.;
		fact_n = 1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				fact_n /= (double)i;
				range = evalrange(x);
				r += fact_n * exp(range) * hn;
			}else {
				hn *= h;
				fact_n /= (double)i;
				r += fact_n * ea * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa sqrt(const psa& x) {
		T a, xn, xn2, sqrt_a, range;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		sqrt_a = sqrt(a);

		r = sqrt_a;
		hn = 1.;
		xn = 1./(2. * sqrt_a);
		range = evalrange(x);
		xn2 = 1./(2. * sqrt(range));
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				r += xn2 * hn;
			} else {
				hn *= h;
				r += xn * hn;
				xn *= (1./2. - i) / a / (i + 1.);
				xn2 *= (1./2. - i) / range / (i + 1.);
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa log(const psa& x) {
		T a, xn, xn2, range;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		r = log(a);
		hn = 1.;
		xn = -1.;
		range = evalrange(x);
		xn2 = -1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				xn2 = -xn2 / range;
				r += xn2 / (double)i * hn;
			} else {
				hn *= h;
				xn = -xn / a;
				xn2 = -xn2 / range;
				r += xn / (double)i * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa sinh (const psa& x) {
		T a, xn, xn2, range, fact_n, table[2], tmp;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		table[0] = sinh(a);
		table[1] = cosh(a);

		r = table[0];
		hn = 1.;
		fact_n = 1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				fact_n /= (double)i;
				range = evalrange(x);
				switch (i%2) {
				case 0: tmp = sinh(range); break;
				case 1: tmp = cosh(range); break;
				}
				r += fact_n * tmp * hn;
			} else {
				hn *= h;
				fact_n /= (double)i;
				r += fact_n * table[i % 2] * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa cosh (const psa& x) {
		T a, xn, xn2, range, fact_n, table[2], tmp;
		psa h, hn, r;
		int i;
		int old_size;
		bool recover_uh = false;
		bool recover_rh = false;

		a = x.v(0);
		// h = x - a;
		h = x; h.v(0) = 0.;

		table[0] = cosh(a);
		table[1] = sinh(a);

		r = table[0];
		hn = 1.;
		fact_n = 1.;
		if (use_history() == true) {
			old_size = history().front().v.size();
		}
		for (i=1; i<x.v.size(); i++) {
			if (use_history() == true && i >= old_size) {
				use_history() = false;
				recover_uh = true;
			}
			if (mode() == 2 && i == x.v.size() - 1) {
				if (record_history() == true) {
					record_history() = false;
					recover_rh = true;
				}
				hn *= h;
				fact_n /= (double)i;
				range = evalrange(x);
				switch (i%2) {
				case 0: tmp = cosh(range); break;
				case 1: tmp = sinh(range); break;
				}
				r += fact_n * tmp * hn;
			} else {
				hn *= h;
				fact_n /= (double)i;
				r += fact_n * table[i % 2] * hn;
			}
		}
		if (recover_uh) use_history() = true;
		if (recover_rh) record_history() = true;

		return r;
	}

	friend psa pow(const psa& x, int y) {
		psa r, xp;
		int a, tmp;

		if (y == 0) return psa(1.);

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

	friend psa pow(const psa& x, const psa& y) {
		return exp(y * log(x));
	}

	friend std::ostream& operator<<(std::ostream& s, const psa& x) {
		int i;
		int n = x.v.size();
		s << '[';
		s << x.v(0);
		for (i=1; i<n; i++) {
			s << ',';
			s << x.v(i);
		}
		s << ']';
		return s;
	}

	friend psa integrate(const psa& x) {
		int i;
		int s = x.v.size();
		psa r;

		r.v.resize(s+1);

		r.v(0) = 0.;
		for (i=1; i<=s; i++) {
			r.v(i) = x.v(i-1) / (double)i;
		}
		return r;
	}

	friend psa setorder(const psa& x, int n) {
		int i;
		int s = x.v.size();
		psa r;

		r.v.resize(n+1);

		if (n+1 >= s) {
			for (i=0; i<s; i++) r.v(i) = x.v(i);
			for (i=s; i<n+1; i++) r.v(i) = 0.;
			return r;
		}

		for (i=0; i<n; i++) r.v(i) = x.v(i);
		if (psa::mode() == 1) {
			r.v(n) = x.v(n);
		} else {
			r.v(n) = polyrange(x.v, n, s-1, psa::domain());
		}

		return r;
	}

	friend T evalrange(const psa& x) {
		int s = x.v.size();

		return polyrange(x.v, 0, s-1, psa::domain());
	}

	friend T eval(const psa& x, const T& a) {
		int s = x.v.size();

		return polyrange(x.v, 0, s-1, a);
	}

	/*
	 *  evaluate { p[x] + p[x+1]t + ... p[y]t^(y-x) | a \in d }
	 */
	static T inline polyrange (const ub::vector<T>& p, int x, int y, const T& d)
	{
		int i;
		T r;

		r = p(y);

		for (i=y-1; i>=x; i--) {
			r = r * d + p(i);
		}

		return r;
	}
};

} // namespace kv

#endif //PSA_HPP
