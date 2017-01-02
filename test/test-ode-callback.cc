#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>
#include <kv/ode-qr.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Lorenz {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};


namespace kv {

template <class T> void callback_sample(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) {
	std::cout << "t: " << end << "\n";
	std::cout << x_e << "\n";
}

} // namespace kv

namespace kv {

template <class T, int N> void callback_sample2(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) {
	int i, j;
	interval<T> t;
	ub::vector< interval<T> > y;
	psa< interval <T> > tmp;
	int s = result.size();
	y.resize(s);
	for (i=1; i<=N; i++) {
		t = (end - start)  * (double)i / (double)N;
		for (j=0; j<s; j++) {
			tmp = result(j);
			y(j) = eval(tmp, t);
		}
		std::cout << "t: " << start + t << "\n";
		std::cout << y << "\n";
	}
}

} // namespace kv


int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	ix = x;
	end = 1.;
	r = kv::odelong_maffine(Lorenz(), ix, itvd(0.), end, kv::ode_param<double>(), kv::callback_sample<double>);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 1.;
	r = kv::odelong_qr(Lorenz(), ix, itvd(0.), end, kv::ode_param<double>(), kv::callback_sample<double>);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 1.;
	r = kv::odelong_maffine(Lorenz(), ix, itvd(0.), end, kv::ode_param<double>(), kv::callback_sample2<double, 3>);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 1.;
	r = kv::odelong_qr(Lorenz(), ix, itvd(0.), end, kv::ode_param<double>(), kv::callback_sample2<double, 3>);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
