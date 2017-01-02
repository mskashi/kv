#include <iostream>
#include <limits>

#include <kv/ode-qr-lohner.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};


int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itv> ix;
	ub::vector< kv::autodif<itv> > dx;
	int r;

	itv end;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	ix = x;
	end = 1.;
	r = kv::odelong_qr_lohner(Lorenz(), ix, itv(0.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	dx = kv::autodif<itv>::init(ix);
	end = 1.;
	r = kv::odelong_qr_lohner(Lorenz(), dx, itv(0.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << dx << "\n";
		std::cout << end << "\n";
	}
}
