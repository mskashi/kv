#include <kv/allsol.hpp>
#include <kv/psa.hpp>
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

	Legendre(int n) : n(n){
	}

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		int i;
		T tmp1, tmp2, tmp3;
		ub::vector<T> y(1);

		if (n == 0) {
			y(0) = 1.;
			return y;
		}
		if (n == 1) {
			y(0) = x(0);
			return y;
		}
		tmp1 = 1.;
		tmp2 = x(0);
		for (i=2; i<=n; i++) {
			tmp3 = ((2*i-1) * x(0) * tmp2 - (i-1) * tmp1) / i;
			tmp1 = tmp2;
			tmp2 = tmp3;
		}
		y(0) = tmp3;
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

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		ub::vector<T> y(1);

		y(0) = eval(p, x(0));

		return y;
	}
};

int main()
{
	ub::vector<itv> I(1);
	std::cout.precision(17);

	I(0) = itv(-1, 1);

	// allsol(Legendre(20), I, 2);
	allsol(Legendre_psa<itv>(20), I, 2);
}
