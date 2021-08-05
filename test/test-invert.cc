#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/random.hpp>
#include <ctime>
#include <chrono>

#include <kv/matrix-inversion.hpp>

int main()
{
	boost::numeric::ublas::matrix<double> a(2, 2);
	boost::numeric::ublas::matrix<double> b;

	std::chrono::system_clock::time_point t;

	a(0, 0) = 1.; a(0, 1) = 2.;
	a(1, 0) = 3.; a(1, 1) = 4.;

	kv::invert(a, b);

	std::cout << b << "\n";

	int i, j, n;

	n = 2000;

	boost::numeric::ublas::matrix<double> c(n, n);
	using namespace boost;
	variate_generator< mt19937, uniform_real<double> > rand (mt19937(time(0)), uniform_real<double>(-1., 1.));

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			c(i, j) = rand();
		}
	}

	t = std::chrono::system_clock::now();
	kv::invert(c, b);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	boost::numeric::ublas::vector<double> d(n);

	t = std::chrono::system_clock::now();
	kv::linear_equation(c, d, d);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	t = std::chrono::system_clock::now();
	kv::mm_mult(c, b, b);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
}
