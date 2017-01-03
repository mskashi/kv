#include <iostream>
#include <limits>
#include <kv/ode-maffine.hpp>

#include "ivp-example.hpp"


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


int main()
{
	ub::vector<itv> x;
	itv start, end;
	kv::ode_param<double> p;
	bool r;

	std::cout.precision(17);


	SimpleHarmonic p01;

	p01.initial_value(x);
	p01.start_time(start);
	p01.stop_time(end);

	std::cout << "SimpleHarmonic\n";
	r = odelong_maffine(p01, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	Lorenz p02;

	p02.initial_value(x);
	p02.start_time(start);
	p02.stop_time(end);

	std::cout << "Lorenz\n";
	r = odelong_maffine(p02, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	VdP p03;

	p03.initial_value(x);
	p03.start_time(start);
	p03.stop_time(end);

	std::cout << "VdP\n";
	r = odelong_maffine(p03, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	Nobi p04;

	p04.initial_value(x);
	p04.start_time(start);
	p04.stop_time(end);

	std::cout << "Nobi\n";
	r = odelong_maffine(p04, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	QuadTest1 p05;

	p05.initial_value(x);
	p05.start_time(start);
	p05.stop_time(end);

	std::cout << "QuadTest1\n";
	r = odelong_maffine(p05, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	QuadTest2 p06;

	p06.initial_value(x);
	p06.start_time(start);
	p06.stop_time(end);

	std::cout << "QuadTest2\n";
	r = odelong_maffine(p06, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
