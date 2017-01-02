#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "dd.hpp"

typedef kv::dd dd;

namespace ub = boost::numeric::ublas;

int main()
{
	ub::matrix<kv::dd> L(2,2), b(2,2);
	ub::permutation_matrix<> pm(2);
	std::cout.precision(40);

	L(0,0) = -10.;
	L(0,1) = 10.;
	L(1,0) = 1.;
	L(1,1) = -1.;

	std::cout << 1/L(0,0) << "\n";
	std::cout << L(1,0)*(1/L(0,0)) << "\n";
	std::cout << L(1,0)*(1/L(0,0))*L(0,1) << "\n";
	std::cout << L(1,1) - L(1,0)*(1/L(0,0))*L(0,1) << "\n";

	if (ub::lu_factorize(L, pm) != 0) {
		std::cout << "singular\n";
	} 
}
