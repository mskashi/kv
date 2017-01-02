#include <iostream>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <kv/lp.hpp>

namespace ub = boost::numeric::ublas;

int main()
{
	double r;
	
	ub::vector<double> objfunc;
	std::list< ub::vector<double> > constraints;
	ub::vector<double> tmp;

	objfunc.resize(3);
	objfunc(0) = 0.;
	objfunc(1) = -2.;
	objfunc(2) = -1.;

	tmp.resize(3);
	// tmp(0) = -5.;
	tmp(0) = -5.2;
	tmp(1) = 1.;
	tmp(2) = -1.;
	constraints.push_back(tmp);

	tmp(0) = -10.;
	tmp(1) = 1.;
	tmp(2) = 1.;
	constraints.push_back(tmp);

	std::cout.precision(17);

	r = kv::lp_minimize_verified(objfunc, constraints, -1);
	std::cout << r << "\n";
	r = kv::lp_minimize(objfunc, constraints);
	std::cout << r << "\n";
	r = kv::lp_minimize_verified(objfunc, constraints, 1);
	std::cout << r << "\n";
}
