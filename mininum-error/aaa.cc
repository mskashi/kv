#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "interval.hpp"
#include "dd.hpp"
#include "rdd.hpp"


namespace ub = boost::numeric::ublas;

int main() {
	kv::interval< kv::dd > hoge;

	ub::matrix< kv::dd > L;

	ub::lu_factorize(L);
}
