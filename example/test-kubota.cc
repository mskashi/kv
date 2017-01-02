#include <iostream>
#include <limits>

#include <kv/ode-maffine2.hpp>


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

	std::cout.precision(17);

	x.resize(4);
	x(0) = -0.9;
	x(1) = -1.4;
	x(2) = 0.01;
	x(3) = 0.4;

	Kepler f;

	ix = x;
	end = std::numeric_limits<double>::infinity();
	r = kv::odelong_maffine2(f, ix, itvd(0.), end, kv::ode_param<double>().set_verbose(1).set_restart_max(10).set_ep_reduce(1000).set_ep_reduce_limit(1100));
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
