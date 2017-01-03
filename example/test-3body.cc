#include <iostream>
#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct ThreeBody {
	template <class T> T pow23(T x, T y) {
		T tmp;
		tmp = x*x + y*y;
		return tmp * sqrt(tmp);
	}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(12);

		T x1 = x(0);
		T y1 = x(1);
		T x2 = x(2);
		T y2 = x(3);
		T x3 = x(4);
		T y3 = x(5);

		y(0) = x(6);
		y(1) = x(7);
		y(2) = x(8);
		y(3) = x(9);
		y(4) = x(10);
		y(5) = x(11);

		T d12 = pow23(x1-x2,y1-y2);
		T d13 = pow23(x1-x3,y1-y3);
		T d23 = pow23(x2-x3,y2-y3);

		y(6) = (-12.*(x1-x2)/d12 - 15.*(x1-x3)/d13) / 3.;
		y(7) = (-12.*(y1-y2)/d12 - 15.*(y1-y3)/d13) / 3.;
		y(8) = (12.*(x1-x2)/d12 - 20.*(x2-x3)/d23) / 4.;
		y(9) = (12.*(y1-y2)/d12 - 20.*(y2-y3)/d23) / 4.;
		y(10) = (15.*(x1-x3)/d13 + 20.*(x2-x3)/d23) / 5.;
		y(11) = (15.*(y1-y3)/d13 + 20.*(y2-y3)/d23) / 5.;

		return y;
	}
};

int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itv> ix;
	bool r;

	std::cout.precision(17);

	itv end;

	x.resize(12);
	x(0) = 1.; x(1) = 3.;
	x(2) = -2.; x(3) = -1.;
	x(4) = 1.; x(5) = -1.;
	x(6) = x(7) = x(8) = x(9) = x(10) = x(11) = 0.;


	ThreeBody f;

	ix = x;
	end = 70.;
	r = kv::odelong_maffine(f, ix, itv(0.), end, kv::ode_param<double>().set_verbose(1).set_restart_max(10));
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
