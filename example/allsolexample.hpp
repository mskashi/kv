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
			x(i).assign(-8., 8.);
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
			x(i).assign(-8., 8.);
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
			x(i).assign(-10., 10.);
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
			x(i).assign(-10., 10.);
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
			x(i).assign(-1e20, 1e20);
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
			x(i).assign(0., 1e10);
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
			x(i).assign(0.01, 100.);
		}
	}
};

// 篠原能材: 数値解析の基礎 問3.8.2
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
			x(i).assign(-5., 5.);
		}
	}
};

// 篠原能材: 数値解析の基礎 問3.8.3
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
			x(i).assign(-3., 3.);
		}
	}
};


// 篠原能材: 数値解析の基礎 演習問題3.8
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
		x(0).assign(-1.5, 2.);
		x(1).assign(-0.6, 3.);
		x(2).assign(-1.5, 2.5);
		x(3).assign(-0.5, 1.9);
		x(4).assign(-1, 1.);
	}
};


/*
  http://nlab.ee.tokushima-u.ac.jp/nishio/Pub-Data/WORK/W153.pdf
  SPICEを用いた複数個の動作点解析 -全解探索への挑戦-
  に出てきた例題。元は、
  http://200.13.98.241/~martin/syop/tareas3/kuno_seader_homotopy.pdf
  Computing All Real Solutions to Systems of Nonlinear Equations
  with a Global Fixed-Point Homotopy
  に出てくる、Himmelblauを少し変形したもの。
  7つ解があるが、ホモトピー法でうまく解が見つからないことがあるらしい。
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
			x(i).assign(-1e8, 1e8);
		}
	}
};

//
// 2013/07/25のゼミで平八郎君が適当に作った例題
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
			x(i).assign(-1000., 1000.);
		}
	}
};


//
// 山村先生の論文によく出てくる区分線形な方程式
//

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
			x(i).assign(-2.5, 2.5);
		}
	}
};
