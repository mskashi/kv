#include <iostream>
#include <limits>

#include "ode-maffine.hpp"
#include "odescale.hpp"


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Kepler {
	public:

	template <class T> T pow23(T x, T y) {
			T tmp;
			tmp = x*x + y*y;
			return tmp * sqrt(tmp);
	}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(4);

		T px = x(0);
		T py = x(1);
		T qx = x(2);
		T qy = x(3);

		T d = pow23(qx,qy);

		y(0) = -qx / d;
		y(1) = -qy / d;
		y(2) = px;
		y(3) = py;

		return y;
	}
};

int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	x.resize(4);
	x(0) = 0.3;
	x(1) = 0.2;
	x(2) = 3.;
	x(3) = 0.;

	std::cout.precision(17);

	Kepler f;
	kv::Odescale<Kepler> g(f);

	ix = x;
	// end = std::numeric_limits<double>::infinity();
	end = 16.;
	r = kv::odelong_maffine(f, ix, itvd(0.), end, 15, 2, 1);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(5);
	for (i=0; i<4; i++) ix(i) = x(i);
	ix(4) = 0.; // start time of original problem
	// end = std::numeric_limits<double>::infinity();
	end = 25.3;
	r = kv::odelong_maffine(g, ix, itvd(0.), end, 15, 2, 1);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
