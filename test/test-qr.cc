#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <ctime>

#include <kv/qr.hpp>

namespace ub = boost::numeric::ublas;

int main()
{
	int i, j;
	ub::matrix<double> a, q, r;

	a.resize(2, 2);

	a(0, 0) = 1.;
	a(0, 1) = 2.;
	a(1, 0) = 3.;
	a(1, 1) = 4.;

	kv::qr(a, q, r);

	std::cout << q << "\n";
	std::cout << r << "\n";
	std::cout << a - prod(q, r) << "\n";
	std::cout << prod(q, trans(q)) << "\n";

	for (i=0; i<10; i++) {
		kv::qr(a, q, r);
		a = prod(r, q);
		std::cout << a << "\n";
	}

	a.resize(5, 5);
	srand(time(NULL));
	for (i=0; i<5; i++) {
		for (j=0; j<5; j++) {
			a(i, j) = (double)rand();
		}
	}

	kv::qr(a, q, r);

	// std::cout << q << "\n";
	// std::cout << r << "\n";
	std::cout << a - prod(q, r) << "\n";
	std::cout << prod(q, trans(q)) << "\n";
}
