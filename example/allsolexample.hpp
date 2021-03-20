#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>

namespace ub = boost::numeric::ublas;

struct Matsu1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = (x(0)-1.) * (x(0)-1.) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-8., 8.);
		}
	}
};

struct Matsu2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1) + 1.;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-8., 8.);
		}
	}
};

struct NoSol {

	double param;

	NoSol(double param = 1e-5) : param(param) {
	}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) - x(1);
		y(1) = x(0) * x(0) - x(1) + param;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-10., 10.);
		}
	}
};

struct BadCond {

	double param;

	BadCond (double param = 1e-7) : param(param) {
	}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) - x(1);
		y(1) = (1. - param) * x(0) * x(0) - x(1) + param;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-10., 10.);
		}
	}
};

struct Hansen1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(1);

		y(0) = x(0)*x(0)*x(0)*x(0) - 12.*x(0)*x(0)*x(0) + 47.*x(0)*x(0) - 60.*x(0);

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(1);
		for (i=0; i<1; i++) {
			x(i) = kv::interval<T>(-1e20, 1e20);
		}
	}
};

struct Burden {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * (4. - 0.0003 * x(0) - 0.0004 * x(1));
		y(1) = x(1) * (2. - 0.0002 * x(0) - 0.0001 * x(1));

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(0., 1e10);
		}
	}
};

struct GE1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = -1. / (1. + x(0)) + x(1) / (1. + x(1));
		y(1) = 1. / (1. + x(0)) - x(1) / (x(0) + x(1));

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(0.01, 100.);
		}
	}
};

// Yoshitane Shinohara: Suuchikaiseki no Kiso, q. 3.8.2
struct Shinohara1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);
		T p, q;

		p = x(0);
		q = x(1);

		y(0) = p*p*p*p*p - 10.*p*p*p*q*q + 5.*p*q*q*q*q
		       - 2.*p*p*p*p + 12.*p*p*q*q - 2.*q*q*q*q
		       + 10.*p*p*p - 30.*p*q*q - 9.*p + 3.;
		y(1) = q*q*q*q*q - 10.*p*p*q*q*q + 5.*p*p*p*p*q
		       - 8.*p*p*p*q + 8.*p*q*q*q + 30.*p*p*q
		       -10.*q*q*q - 9.*q;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-5., 5.);
		}
	}
};

// Yoshitane Shinohara: Suuchikaiseki no Kiso, q. 3.8.3
struct Shinohara2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);
		T p, q;

		p = x(0);
		q = x(1);

		y(0) = p*p*p*p*p - 10.*p*p*p*q*q + 5.*p*q*q*q*q
		       - 3.*p*p*p*p + 18.*p*p*q*q - 3.*q*q*q*q
		       - 2.*p*p*p + 6.*p*q*q + 3.*p*p*q - q*q*q
		       + 12.*p*p - 12.*q*q - 10.*p*q - 8.*p + 8.*q;
		y(1) = 5.*p*p*p*p*q - 10.*p*p*q*q*q + q*q*q*q*q
		       - 12.*p*p*p*q + 12.*p*q*q*q - p*p*p + 3.*p*q*q
		       - 6.*p*p*q + 2.*q*q*q + 5.*p*p - 5.*q*q
		       + 24.*p*q - 8.*p - 8.*q + 4.;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-3., 3.);
		}
	}
};


// Yoshitane Shinohara: Suuchikaiseki no Kiso, ex. 3.8
struct Shinohara3 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(5);
		T p, q, r, s, t;

		p = x(0);
		q = x(1);
		r = x(2);
		s = x(3);
		t = x(4);

		y(0) = p*p*p - 2.*p*q + r + 0.75*p + 1.;
		y(1) = p*p*q - q*q - p*r +s + 0.75*q + 0.25;
		y(2) = p*p*r - p*s - q*r + t + 0.75*r + 0.75;
		y(3) = p*p*s - p*t - q*s + 0.75*s;
		y(4) = p*p*t - q*t + 0.75*t - 0.25;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(5);
		x(0) = kv::interval<T>(-1.5, 2.);
		x(1) = kv::interval<T>(-0.6, 3.);
		x(2) = kv::interval<T>(-1.5, 2.5);
		x(3) = kv::interval<T>(-0.5, 1.9);
		x(4) = kv::interval<T>(-1, 1.);
	}
};


/*
  Problem taken from
  http://nlab.ee.tokushima-u.ac.jp/nishio/Pub-Data/WORK/W153.pdf, 
  http://200.13.98.241/~martin/syop/tareas3/kuno_seader_homotopy.pdf
  Computing All Real Solutions to Systems of Nonlinear Equations
  with a Global Fixed-Point Homotopy
  This equation has 7 solutions and difficult to find all solution
  by homotopy method.
 */

struct ModifiedHimmelblau {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = 2.*x(0)*x(0)*x(0) + 2.*x(0)*x(1) - 22.*x(0) + x(1)*x(1) + 13.;
		y(1) = x(0)*x(0) + 2.*x(0)*x(1) + 2.*x(1)*x(1)*x(1) - 14.*x(1) + 9.;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-1e8, 1e8);
		}
	}
};

//
// made by Heihachiro Yoshii on 2013/07/25
//

struct Heihachiro {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = 2. * x(0) * x(0) * x(1) - 1.;
		y(1) = x(0) + 0.5 * x(1) * x(1) - 2.;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i) = kv::interval<T>(-1000., 1000.);
		}
	}
};


struct Yamamura2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int n = x.size();
		ub::vector<T> y(n);
		T s;
		int i;

		s = 0.;
		for (i=0; i<n; i++) {
			s += pow(x(i), 3);
		}

		for (i=0; i<n; i++) {
			y(i) = x(i) - (s + (i + 1.)) / (2. * n);
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i) = kv::interval<T>(-2.5, 2.5);
		}
	}
};


//
// K. Meintjes and A. P. Morgan:
// Chemical Equilibrium Systems as Numerical Test Problmes,
// ACM Transactions on Mathematical Software, 16(2):143, 1990.
//

struct HydroCarbon {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(5);
		// phisical constants
		static T K5 = kv::constants<T>::str("1.930e-1");
		static T K6 = kv::constants<T>::str("2.597e-3");
		static T K7 = kv::constants<T>::str("3.448e-3");
		static T K8 = kv::constants<T>::str("1.799e-5");
		static T K9 = kv::constants<T>::str("2.155e-4");
		static T K10 = kv::constants<T>::str("3.846e-5");
		// parameters
		static T R = T(10);
		static T p = T(40);
		// constants
		static T R5 = K5;
		static T R6 = K6 * pow(p, -0.5);
		static T R7 = K7 * pow(p, -0.5);
		static T R8 = K8 * pow(p, -1);
		static T R9 = K9 * pow(p, -0.5);
		static T R10 = K10 * pow(p, -1);

		y(0) = x(0) * x(1) + x(0) - 3 * x(4);
		y(1) = 2 * x(0) * x(1) + x(0) + 2 * R10 * pow(x(1), 2) + x(1) * pow(x(2), 2) + R7 * x(1) * x(2) + R9 * x(1) * x(3) + R8 * x(1) - R * x(4);
		y(2) = 2 * x(1) * pow(x(2), 2) + R7 * x(1) * x(2) + 2 * R5 * pow(x(2), 2) + R6 * x(2) - 8 * x(4);
		y(3) = R9 * x(1) * x(3) + 2 * pow(x(3), 2) - 4 * R * x(4);
		y(4) = x(0) * x(1) + x(0) + R10 * pow(x(1), 2) + x(1) * pow(x(2), 2) + R7 * x(1) * x(2) + R9 * x(1) * x(3) + R8 * x(1) + R5 * pow(x(2), 2) + R6 * x(2) + pow(x(3), 2) - 1;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		x.resize(5);
		int i;
		for (i=0; i<5; i++) {
			x(i) = kv::interval<T>(-100, 100);
		}
	}
};
