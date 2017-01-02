#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>


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
	ub::vector< kv::affine<double> > ax;
	kv::affine<double> dummy;
	int r;

	itv end;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	ax = x;
	end = std::numeric_limits<double>::infinity();
	r = kv::ode_maffine(Lorenz(), ax, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		for (i=0; i<ax.size(); i++) {
			std::cout << to_interval(ax(i)) << "\n";
		}
		std::cout << end << "\n";
	}

	kv::affine<double>::maxnum() = 0;
	ax = x;
	end = 1.;
	r = kv::odelong_maffine(Lorenz(), ax, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		for (i=0; i<ax.size(); i++) {
			std::cout << to_interval(ax(i)) << "\n";
		}
		std::cout << end << "\n";
	}

	ix = x;
	end = 1.;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	dx = kv::autodif<itv>::init(ix);
	end = 1.;
	r = kv::odelong_maffine(Lorenz(), dx, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << dx << "\n";
		std::cout << end << "\n";
	}
}
