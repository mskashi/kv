#include <iostream>
#include <limits>

#include <kv/ode-affine.hpp>
#include <kv/ode-affine-wrapper.hpp>


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
	ub::vector< kv::affine<double> > x;
	int r;
	itv end;
	int i;

	std::cout.precision(17);

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;

	r = kv::odelong_affine(Lorenz(), x, (itv)0., end);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		for (i=0; i<x.size(); i++) {
			std::cout << to_interval(x) << "\n";
		}
		std::cout << end << "\n";
	}


	// reset maxnum of affine to speed up
	kv::affine<double>::maxnum() = 0;

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;

	r = kv::odelong_wrapper(Lorenz(), x, (itv)0., end);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		for (i=0; i<x.size(); i++) {
			std::cout << to_interval(x) << "\n";
		}
		std::cout << end << "\n";
	}
}
