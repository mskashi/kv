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
		ub::vector<T> y(12);

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

	double e = 0.5;

	itvd end;

	x.resize(4);
	x(0) = 0.;
	x(1) = sqrt((1+e)/(1-e));
	x(2) = 1-e;
	x(3) = 0.;

	std::cout.precision(17);

	Kepler f;

#if 0
	ix.resize(13);
	for (i=0; i<12; i++) ix(i) = x(i);
	ix(12) = 0.;
	end = std::numeric_limits<double>::infinity();

	r = kv::odelong_maffine(g, ix, itvd(0.), end, 24, 2, 1);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
#endif

	ix = x;
	end = std::numeric_limits<double>::infinity();
	r = kv::odelong_maffine(f, ix, itvd(0.), end, 15, 2, 1);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
