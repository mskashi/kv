#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>
#include <kv/psa-plot.hpp>


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


namespace kv {

template <class T> struct ode_callback_psaplot : ode_callback<T> {
	matplotlib g;
	int n;

	ode_callback_psaplot(matplotlib g, int n) : g(g), n(n) {}

	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		int i, j;
		int s = result.size();

		for (j=0; j<s; j++) {
			psa< interval<T> >::domain() = interval<T>(0., (end - start).upper());
			psa_plot(result(j), mid(start), g, n);
		}
		return true;
	}
};

} // namespace kv


int main()
{
	int i;
	ub::vector<itv> ix;
	int r;
	kv::matplotlib g;

	itv end;

	ix.resize(3);
	ix(0) = 15; ix(1) = 15.; ix(2) = 36.;
	// ix(0) = itv(15., 15.1); ix(1) = itv(15., 15.1); ix(2) = itv(36., 36.1);

	std::cout.precision(17);

	g.open();
	g.screen(0, -50, 1, 50);

	end = 1.;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_psaplot<double>(g, 10));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	getchar();
}
