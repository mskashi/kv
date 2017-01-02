#include <iostream>
#include <limits>
#include <kv/ode.hpp>

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

int main()
{
	ub::vector<itvd> x;
	int r;
	itvd end;

	std::cout.precision(17);

	x.resize(3);

	//
	// ode(f, x, start, end, order, autostep=true, iter_max=2, result_psa=NULL)
	// 
	// integrate f from start (0) to end (0.1)
	// x : initial vector
	// end time should be given by "writable" interval variable
	// order (24): order of Taylor expansion.
	// autostep: use automatic step size control or not
	// iter_max: maximum number of iterative refinement
	// the result is overwritten to x.
	// return value: 0 - can't calculate verified solution.
	//               1 - verified solution is obtained but end time
	//                   is less than given end time.
	//                   the reference valiable end is overwritten.
	//               2 - verified solution is obtained up to given end time
	// if result_psa != NULL then the verified solution on [start, end] is written
	// to the variable result_psa.
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 0.1;
	r = kv::ode(Lorenz(), x, itvd(0.), end, 24, false);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// too long step size
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 10.;
	r = kv::ode(Lorenz(), x, itvd(0.), end, 24, false);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// enable automatic step size control
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 10.;
	r = kv::ode(Lorenz(), x, itvd(0.), end, 24, true);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// if you want to progress as long as possible, set end to infinity.
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = std::numeric_limits<double>::infinity();
	r = kv::ode(Lorenz(), x, itvd(0.), end, 24, true);
	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// odelong : connect multiple ode (autostep=true)
	// odelong(f, x, start, end, order, iter_max=2, verbose=0)
	// if verbose=1 then become verbose
	// NOTICE: odelong don't use any idea to supress so-called wrapping effect.
	//         so odelong can't calculate long time.
	//         for long time calculation, you should use odelong_maffine or others.
	end = 1.;
	r = kv::odelong(Lorenz(), x, itvd(0.), end, 24);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
