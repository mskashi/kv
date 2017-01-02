#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class DoublePendulum {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(4);
		ub::matrix<T> tm(2,2);
		ub::vector<T> tv(2);
		T a, b, c, d;

		a = 2.;
		b = cos(x(0) - x(1));
		c = b;
		d = 1.;

		tm(0,0) = d;
		tm(0,1) = -b;
		tm(1,0) = -c;
		tm(1,1) = a;
		tm /= (a*d - b*c);

		tv(0) = -9.81 * 2 * sin(x(0)) - sin(x(0) - x(1)) * x(3) * x(3);
		tv(1) = -9.81 * sin(x(1)) + sin(x(0) - x(1)) * x(2) * x(2);

		tv = prod(tm, tv);

		y(0) = x(2);
		y(1) = x(3);
		y(2) = tv(0);
		y(3) = tv(1);

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	itvd p = kv::constants<itvd>::pi();

	ix.resize(4);
	ix(0) = itvd(0.99, 1.01) * 3 * p / 4;
	ix(1) = -11 * p / 20;
	ix(2) = 0.43;
	ix(3) = 0.67;

	ix(1) += ix(0);
	ix(3) += ix(2);

	std::cout.precision(17);

	// end = std::numeric_limits<double>::infinity();
	end = 0.5;
	r = kv::odelong_maffine(DoublePendulum(), ix, itvd(0.), end, kv::ode_param<double>().set_verbose(1));
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		for (i=0; i<ix.size(); i++) {
			std::cout << ix(i) << "\n";
		}
		std::cout << end << "\n";
	}
}
