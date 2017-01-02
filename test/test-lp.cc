#include <iostream>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <kv/lp.hpp>

#ifdef TEST_DD
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
typedef kv::dd fl;
#else
typedef double fl;
#endif


namespace ub = boost::numeric::ublas;

int main()
{
	fl r;
	kv::interval<fl> ri;
	
	ub::vector<fl> objfunc;
	std::list< ub::vector<fl> > constraints;
	ub::vector<fl> tmp;

#ifdef TEST_DD
	std::cout.precision(33);
#else
	std::cout.precision(17);
#endif

	// simple problem

	objfunc.resize(3);
	tmp.resize(3);

	objfunc(0) = 0.;
	objfunc(1) = -2.;
	objfunc(2) = -1.;

	tmp(0) = -5.;
	tmp(1) = 1.;
	tmp(2) = -1.;
	constraints.push_back(tmp);

	tmp(0) = -10.;
	tmp(1) = 1.;
	tmp(2) = 2.;
	constraints.push_back(tmp);

	r = kv::lp_minimize(objfunc, constraints);
	std::cout << r << "\n";

	ri.lower() = kv::lp_minimize_verified(objfunc, constraints, -1);
	ri.upper() = kv::lp_minimize_verified(objfunc, constraints, 1);
	std::cout << ri << "\n";

	// simple problem by akky

	constraints.clear();
	objfunc.resize(4);
	tmp.resize(4);

	objfunc(0) = 0.;
	objfunc(1) = 1.;
	objfunc(2) = 0.;
	objfunc(3) = -1.;

	tmp(0) = -2.;
	tmp(1) = 1.;
	tmp(2) = 0.;
	tmp(3) = 0.;
	constraints.push_back(tmp);

	tmp(0) = 0.;
	tmp(1) = -2.;
	tmp(2) = 1.;
	tmp(3) = 0.;
	constraints.push_back(tmp);

	tmp(0) = -4.;
	tmp(1) = 4.;
	tmp(2) = -1.;
	tmp(3) = 0.;
	constraints.push_back(tmp);

	tmp(0) = 0.;
	tmp(1) = -2.;
	tmp(2) = 0.;
	tmp(3) = 1.;
	constraints.push_back(tmp);

	tmp(0) = 0.;
	tmp(1) = -4.;
	tmp(2) = 2.;
	tmp(3) = 1.;
	constraints.push_back(tmp);

	tmp(0) = -4.;
	tmp(1) = 6.;
	tmp(2) = -2.;
	tmp(3) = -1.;
	constraints.push_back(tmp);

	r = kv::lp_minimize(objfunc, constraints);
	std::cout << r << "\n";

	ri.lower() = kv::lp_minimize_verified(objfunc, constraints, -1);
	ri.upper() = kv::lp_minimize_verified(objfunc, constraints, 1);
	std::cout << ri << "\n";
}
