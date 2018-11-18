#include <kv/allsol.hpp>
#include <kv/psa.hpp>
#include <kv/dka.hpp>
#include <kv/constants.hpp>
#include <kv/kraw-approx.hpp>
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


// calculate all zeros of Legendre polynomial
// using approximate initial value 
// (x_i = cos((i+0.75)/(n+0.5)*pi))   (i=0...n-1)

template <class T, class F>
bool
legendle_zeros(F f, int n, ub::vector< kv::interval<T> > &result) {

	int i, j;
	bool r;
	T app;

	result.resize(n);

	for (i=0; i<n; i++) {
		app = cos((i+T(0.75))/(n+T(0.5)) * kv::constants<T>::pi());
		r = kv::krawczyk_approx(f, app, result(i), 5, 0);
		if (!r) {
			// std::cout << "result(" << i << ") failed\n";
			return false;
		}
	}
	// check all solutions are separated from each other
	for (i=0; i<n-1; i++) {
		for (j=i+1; j<n; j++) {
			if (overlap(result(i), result(j))) {
				return false;
			}
		}
	}
	return true;
}


int main()
{
	std::cout.precision(17);

	int n = 20;

	// function object of Legendre Polynomial using recursive evaluation
	Legendre f(n);

	// function object of Legendre Polynomial using pre-calculation of coefficients
	Legendre_psa<itv> g(n);


	// calculate all zeros of Legendre polynomial 
	// using branch-and-bound algorithm

	std::cout << "[allsol, recursive]\n";
	allsol(f, itv(-1, 1), 2);

	std::cout << "[allsol, pre-calculation]\n";
	allsol(g, itv(-1, 1), 2);

	// calculate all zeros of Legendre polynomial
	// using DKA method and Smith's theorem

	std::cout << "[DKA]\n";

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

	// calculate all zeros of Legendre polynomial
	// using approximate initial value

	ub::vector<itv> result;
	bool r;

	std::cout << "[approximate initial value, recursive]\n";
	r = legendle_zeros(f, n, result);
	if (r) {
		std::cout << result << "\n";
	}

	std::cout << "[approximate initial value, pre-calculation]\n";
	r = legendle_zeros(g, n, result);
	if (r) {
		std::cout << result << "\n";
	}
}
