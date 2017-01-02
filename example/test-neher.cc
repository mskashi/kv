#include <iostream>
#include <limits>

#include <kv/ode-affine.hpp>
#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>


// Neher, Jackson, Nedialkov:
// On Taylor model based integration of ODEs
// SIAM J. Numer. Anal., 45:236-262, 2007.


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;

class Neher {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T>& x, T& t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = x(0) * x(0);

		return y;
	}
};

class Neher2 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T>& x, T& t){
		ub::vector<T> y(3);

		ub::matrix<T> B(3,3);

#if 0
		B(0,0) = -0.4375; B(0,1) = 0.0625; B(0,2) = -0.2651650429;
		B(1,0) = 0.0625; B(1,1) = -0.4375; B(1,2) = -0.2651650429;
		B(2,0) = -0.2651650429; B(2,1) = -0.2651650429; B(2,2) = -0.375;
#endif

		B(0,0) = 0.; B(0,1) = -0.7071067810; B(0,2) = -0.5;
		B(1,0) = 0.7071067810; B(1,1) = 0.; B(1,2) = 0.5;
		B(2,0) = 0.5; B(2,1) = -0.5; B(2,2) = 0.;

#if 0
		B(0,0) = -0.125; B(0,1) = -0.8321067810; B(0,2) = -0.3232233048;
		B(1,0) = 0.5821067810; B(1,1) = -0.125; B(1,2) = 0.6767766952;
		B(2,0) = 0.6767766952; B(2,1) = -0.3232233048; B(2,2) = -0.25;
#endif

		y = prod(B, x);

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	ix.resize(2);
	ix(0) = itvd(0.95, 1.05);
	ix(1) = itvd(-1.05, -0.95);

	std::cout.precision(17);

	end = 2.8;
	r = kv::odelong_maffine2(Neher(), ix, itvd(0.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(3);
	ix(0) = itvd(0.999, 1.001);
	ix(1) = itvd(0.999, 1.001);
	ix(2) = itvd(0.999, 1.001);

	end = 100.;
	r = kv::odelong_maffine2(Neher2(), ix, itvd(0.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
