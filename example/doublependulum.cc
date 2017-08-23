#include <iostream>
#include <limits>
#include <kv/ode-maffine.hpp>
#include <kv/constants.hpp>
#include <kv/ode-callback.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct DoublePendulum {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(4);
		ub::matrix<T> tm(2,2);
		ub::vector<T> tv(2);

		T a, b, c, d, t1, t2;

		static const T m1 = T(1.);
		static const T m2 = T(1.);
		static const T l1 = T(1.);
		static const T l2 = T(1.);
		static const T g = kv::constants<T>::str("9.8");

		t1 = cos(x(0) - x(1));
		t2 = sin(x(0) - x(1));

		a = (m1 + m2) * l1;
		b = m2 * l2 * t1;
		c = l1 * l2 * t1;
		d = l2 * l2;

		tm(0,0) = d;
		tm(0,1) = -b;
		tm(1,0) = -c;
		tm(1,1) = a;
		tm /= (a*d - b*c);

		tv(0) = -m2 * l2 * x(3) * x(3) * t2 - (m1 + m2) * g * sin(x(0));
		tv(1) = l1 * l2 * x(2) * x(2) * t2 - g * l2 * sin(x(1));

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
	ub::vector<itv> ix;
	itv end;
	std::cout.precision(17);

	ix.resize(4);
	ix(0) = 0.75 * kv::constants<itv>::pi();
	ix(1) = 0.75 * kv::constants<itv>::pi();
	ix(2) = 0.;
	ix(3) = 0.;

	// end = std::numeric_limits<double>::infinity();
	end = 15.;
	kv::odelong_maffine(DoublePendulum(), ix, itv(0.), end, kv::ode_param<itv::base_type>(), kv::ode_callback_dense_print<itv::base_type>(itv(0.), itv(pow(2., -7))));
}
