#include <iostream>
#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>

namespace ub = boost::numeric::ublas;


/*
  Test Set for IVP Solvers
   https://www.dm.uniba.it/~testset/testsetivpsolvers/
  Robertson problem
   http://www.dm.uniba.it/~testset/problems/rober.php

  - stiff ODE
  - example such that maffine2 is much faster than maffine
*/

struct Rober {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);
		static T c1 = kv::constants<T>::str("0.04");
		static T c2 = kv::constants<T>::str("1e4");
		static T c3 = kv::constants<T>::str("3e7");

		T t1 = c1 * x(0);
		T t2 = c2 * x(1) * x(2);
		T t3 = c3 * pow(x(1), 2);

		y(0) = -t1 + t2;
		y(1) = t1 - t3 - t2;
		y(2) = t3;

		return y;
	}
};

typedef kv::interval<double> itv;

int main()
{
	int i;
	ub::vector<itv> ix;
	int r;
	kv::ode_param<double> p;

	itv end;

	ix.resize(3);

	/*
	// used by Prof. Bunger
	ix(0) = itv(0.99995, 1.00005);
	ix(1) = itv(0, 0.00000001);
	ix(2) = itv(0, 0.00000001);
	*/
	// standard initial value
	ix(0) = 1.;
	ix(1) = 0.;
	ix(2) = 0.;

	std::cout.precision(17);

	end = 1.;
	// odelong_maffine: fast
	// odelong_maffine2: very fast but cannot calculate derivative w.r.t initial value
	r = kv::odelong_maffine(Rober(), ix, itv(0.), end, p.set_verbose(0).set_restart_max(10).set_order(24));
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
