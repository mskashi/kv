#include <kv/allsol.hpp>
#include <kv/psa.hpp>
#include <kv/dka.hpp>
// #include <kv/dd.hpp>
// #include <kv/rdd.hpp>
// #include <kv/mpfr.hpp>
// #include <kv/rmpfr.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;
// typedef kv::interval<kv::dd> itv;
// typedef kv::interval< kv::mpfr<106> > itv;


// Legendre polynomial using recursion

struct Legendre {
	int n;

	Legendre(int n) : n(n) {}

	template <class T> T operator()(const T& x) {
		int i;
		T tmp1, tmp2, y;

		if (n == 0) return T(1);
		if (n == 1) return x;

		tmp1 = 1;
		tmp2 = x;
		for (i=2; i<=n; i++) {
			y = ((2*i-1) * x * tmp2 - (i-1) * tmp1) / i;
			tmp1 = tmp2;
			tmp2 = y;
		}
		return y;
	}
};


// calculate coefficients of Legendre polynomial in advance.

template <class TT> struct Legendre_psa {
	int n;
	kv::psa<TT> p;

	Legendre_psa(int n) : n(n){
		kv::psa<TT> tmp1, tmp2, tmp3, p1, pt;
		int i;

		p1.v.resize(n+1);
		for (i=0; i<n+1; i++) p1.v(i) = 0;
		p1.v(0) = 1.;

		if (n == 0) {
			p = p1;
			return;
		}

		pt = p1;
		pt.v(0) = 0.;
		pt.v(1) = 1.;

		if (n == 1) {
			p = pt;
			return;
		}

		tmp1 = p1;
		tmp2 = pt;
		for (i=2; i<=n; i++) {
			tmp3 = ((2*i-1) * pt * tmp2 - (i-1) * tmp1) / i;
			tmp1 = tmp2;
			tmp2 = tmp3;
		}

		p = tmp3;
	}

	template <class T> T operator()(const T& x) {
		return eval(p, x);
	}
};

int main()
{
	std::cout.precision(17);

	// calculate zeros of Legendre polynomial
	// using recursive evaluation

	Legendre f(20);
	allsol(f, itv(-1, 1), 2);

	// calculate zeros of Legendre polynomial
	// using pre-calculation

	Legendre_psa<itv> g(20);
	allsol(g, itv(-1, 1), 2);

	// calculate zeros of Legendre polynomial
	// using pre-calculation and DKA method

	ub::vector< kv::complex<itv> > va, vr;
	int s = g.p.v.size();
	va.resize(s);
	for (int i=0; i<s; i++) {
		va(i) = g.p.v(i);
	}
	kv::vdka(va, vr);
	for (int i=0; i<vr.size(); i++) {
        	std::cout << vr(i).real() << "\n";
	}
}
