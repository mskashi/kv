#include <boost/timer.hpp>
#include <kv/allsol.hpp>
#ifdef TEST_DD
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#endif

namespace ub = boost::numeric::ublas;

#ifdef TEST_DD
typedef kv::interval<kv::dd> itv;
#else
typedef kv::interval<double> itv;
#endif


//
// 2-transistor circuit equation which has 5 solutions.
//
// Yusuke Nakaya, Tetsuo Nishi, Shin'ichi Oishi, and Martin Claus
// Numerical Existence Proof of Five Solutions for Certain Two-Transistor Circuit Equations
// Japan J. Indust. Appl. Math. Volume 26, Number 2-3 (2009), 327-336. 
//

struct Nishi {
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		ub::vector<T> y(4);
		int i, j;
		ub::vector<T> f(4);
		ub::matrix<T> mT(4,4);
		ub::matrix<T> G(4,4);
		ub::vector<T> J(4);
		T rb, rc, vcc, af, ar, is, vt, vs, gb, gc;

		rb = 10000.;
		rc = 5000.;
		vcc = -5.;
		af = 0.99;
		ar = 0.5;
#ifndef CASEB
		// case (a)
		is = 1e-9;
		vt = 0.053;
		vs = -0.64;
#else
		// case (b)
		is = 1e-6;
		vt = 0.102;
		vs = -0.44;
#endif

		gb = 1. / rb;
		gc = 1. / rc;

		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				mT(i,j) = 0.;
			}
		}
		mT(0,0) = 1.; mT(0,1) = -ar;
		mT(1,0) = -af; mT(1,1) = 1.;
		mT(2,2) = 1.; mT(2,3) = -ar;
		mT(3,2) = -af, mT(3,3) = 1.;

		G(0,0) = 2. * gb + gc;
		G(0,1) = -(gb + gc);
		G(0,2) = -2. * gb;
		G(0,3) = gb;
		G(1,0) = -(gb + gc);
		G(1,1) = gb + gc;
		G(1,2) = gb;
		G(1,3) = 0.;
		G(2,0) = -2. * gb;
		G(2,1) = gb;
		G(2,2) = 2. * gb + gc;
		G(2,3) = -(gb + gc);
		G(3,0) = gb;
		G(3,1) = 0.;
		G(3,2) = -(gb + gc);
		G(3,3) = gb + gc;

		J(0) = gc * vcc;
		J(1) = gb * vs - gc * vcc;
		J(2) = gc * vcc;
		J(3) = gb * vs - gc * vcc;

		f(0) = is / af * (exp(x(0) / vt) - 1.);
		f(1) = is / ar * (exp(x(1) / vt) - 1.);
		f(2) = is / af * (exp(x(2) / vt) - 1.);
		f(3) = is / ar * (exp(x(3) / vt) - 1.);

		y = prod(mT, f) + prod(G, x) + J;

		return y;
	}
};


int main()
{
	int i;
	boost::timer t;
	ub::vector<itv> I;

	#ifdef TEST_DD
	std::cout.precision(34);
	#else
	std::cout.precision(17);
	#endif

	I.resize(4);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	kv::allsol(Nishi(), I);
	std::cout << t.elapsed() << " sec\n";
}
