#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>
#include <kv/ode-affine.hpp>
#include <kv/ode-qr.hpp>
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


namespace kv {

template <class T> struct ode_callback_sample : ode_callback<T> {
	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		std::cout << "t: " << end << "\n";
		std::cout << x_e << "\n";
		return true;
	}
};

} // namespace kv

namespace kv {

template <class T> struct ode_callback_sample2 : ode_callback<T> {
	int n;

	ode_callback_sample2(int n) :n(n) {}

	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		int i, j;
		interval<T> t;
		ub::vector< interval<T> > y;
		psa< interval <T> > tmp;
		int s = result.size();
		y.resize(s);
		for (i=1; i<=n; i++) {
			t = (end - start)  * (double)i / (double)n;
			for (j=0; j<s; j++) {
				tmp = result(j);
				y(j) = eval(tmp, t);
			}
			std::cout << "t: " << start + t << "\n";
			std::cout << y << "\n";
		}
		return true;
	}
};

} // namespace kv


int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itv> ix;
	int r;

	itv end;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	// using ode_callback_sample

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample<double>());
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine2(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample<double>());
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 0.125;
	r = kv::odelong_affine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample<double>());
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 0.125;
	r = kv::odelong_qr(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample<double>());
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix = x;
	end = 0.125;
	r = kv::odelong_qr_lohner(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample<double>());
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	// using ode_callback_sample2

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_sample2<double>(3));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	// test for ode_callback_dense_print in ode-callback.hpp

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_dense_print<double>(itv(0.), itv(0.125)));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	std::list< kv::interval<double> > start_list;
	std::list< kv::interval<double> > end_list;
	std::list< ub::vector< kv::psa< kv::interval<double> > > > psa_list;

	// test for ode_callback_list in ode-callback.hpp

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_list<double>(start_list, end_list, psa_list));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	std::list< kv::interval<double> >::iterator ps, pe;
	std::list< ub::vector< kv::psa< kv::interval<double> > > >::iterator pp;
	ps = start_list.begin();
	pe = end_list.begin();
	pp = psa_list.begin();

	while (ps != start_list.end()) {
		std::cout << *ps << "\n";
		std::cout << *pe << "\n";
		std::cout << *pp << "\n";
		ps++; pe++; pp++;
	}

	// test for ode_callback_dense_list in ode-callback.hpp

	std::list< kv::interval<double> > time_list;
	std::list< ub::vector< kv::interval<double> > > value_list;

	ix = x;
	end = 0.125;
	r = kv::odelong_maffine(Lorenz(), ix, itv(0.), end, kv::ode_param<double>(), kv::ode_callback_dense_list<double>(itv(0.), itv(0.125), time_list, value_list));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	std::list< kv::interval<double> >::iterator pt;
	std::list< ub::vector< kv::interval<double> > >::iterator pv;
	pt = time_list.begin();
	pv = value_list.begin();

	while (pt != time_list.end()) {
		std::cout << *pt << "\n";
		std::cout << *pv << "\n";
		pt++; pv++;
	}
}
