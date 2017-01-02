#include <iostream>
#include <limits>
#include <kv/ode.hpp>

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
	ub::vector<itv> x;
	itv end;
	kv::ode_param<double> p; // parameter for ode solver
	int r;

	std::cout.precision(17);

	x.resize(3);

	//
	// ode(f, x, start, end, ode_param, result_psa=NULL)
	// 
	// integrate f from start to end
	// x : initial vector
	// "x" and "end" must be given by *writable* variable because
	// the result is returned by *overwrite* "x" and "end".
	//
	// effective parameters:
	//   order: order of Taylor expansion.
	//   autostep: use automatic step size control or not
	//   epsilon: relative error tolerance of each step
	//   iteration: maximum number of iterative refinement
	//   restart_max: maximum number of restart when verification failed
	// 
	// return value: 0 - can't calculate verified solution.
	//               1 - verified solution is obtained but end time
	//                   is less than given end time.
	//                   the valiable end is overwritten.
	//               2 - verified solution is obtained up to given end time.
	// if result_psa != NULL then the verified solution on [start, end] is written to the variable result_psa.


	// automatic step size control (default)

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;

	r = kv::ode(Lorenz(), x, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// fixed step size

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 0.1;
	p = kv::ode_param<double>().set_autostep(false);

	r = kv::ode(Lorenz(), x, itv(0.), end, p);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// too long step size

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 10.;
	p = kv::ode_param<double>().set_autostep(false);

	r = kv::ode(Lorenz(), x, itv(0.), end, p);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// if you want to progress as long as possible, set end to infinity.

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = std::numeric_limits<double>::infinity();

	r = kv::ode(Lorenz(), x, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// set many parameters

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;
	p = kv::ode_param<double>().set_order(15).set_iteration(3).set_epsilon(1e-8);

	r = kv::ode(Lorenz(), x, itv(0.), end, p);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	// odelong : connect multiple ode (autostep=true)
	// odelong(f, x, start, end, order, ode_param)
	//
	// effective parameters:
	//   order: order of Taylor expansion.
	//   epsilon: used for deciding error tolerance of each step
	//   iteration: maximum number of iterative refinement
	//   restart_max: maximum number of restart when verification failed
	//   verbose: if verbose==1 then become verbose
	//
	// NOTICE: odelong don't use any idea to supress so-called wrapping effect.
	//         so odelong can't calculate long time.
	//         for long time calculation, you should use odelong_maffine or others.
	end = 1.;
	r = kv::odelong(Lorenz(), x, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
