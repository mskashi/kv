#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(1);
		y(0) = x(0) * (1 - x(0));
		return y;
	}
};

namespace kv {
template <class T> struct ode_callback_stop : ode_callback<T> {
	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		if (x_e(0) >= 1 - std::pow(2., -5)) return false;
		return true;
	}
};
}

int main()
{
	int i;
	ub::vector<itv> ix;
	int r;
	itv end;

	std::cout.precision(17);

	ix.resize(1);
	ix(0) = std::pow(2., -5);
	end = std::numeric_limits<double>::infinity();

	r = kv::odelong_maffine(Func(), ix, itv(0.), end, kv::ode_param<double>().set_verbose(1), kv::ode_callback_stop<double>());
	if (r == 3) {
		std::cout << end << "\n";
		std::cout << ix << "\n";
	}
}
