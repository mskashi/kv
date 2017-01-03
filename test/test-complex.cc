// test program for complex.hpp

#include <kv/complex.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#ifdef TEST_D
typedef kv::complex<double> cp;
#else
 #ifdef TEST_DD
typedef kv::complex<kv::dd> cp;
 #else
  #ifdef TEST_ID
typedef kv::complex< kv::interval<double> > cp;
  #else
typedef kv::complex< kv::interval<kv::dd> > cp;
  #endif // TEST_ID
 #endif // TEST_DD
#endif // TEST_D

int main()
{
	#if defined(TEST_D) || defined(TEST_ID)
	std::cout.precision(17);
	#else
	std::cout.precision(33);
	#endif

	cp x, y, z;

	x = cp(1., 2.);
	y = cp(3., 4.);
	// below is possible because int->double conversion allowed
	x = kv::complex<int>(1, 2);
	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";
	std::cout << cp::i() << "\n";
	std::cout << abs(x) << "\n";
	std::cout << arg(x) << "\n";
	std::cout << conj(x) << "\n";
	std::cout << sqrt(x) << "\n";
	std::cout << pow(x, y) << "\n";
	std::cout << pow(2., y) << "\n";
	std::cout << pow(x, 2.) << "\n";
	std::cout << exp(x) << "\n";
	std::cout << log(x) << "\n";
	std::cout << sin(x) << "\n";
	std::cout << cos(x) << "\n";
	std::cout << tan(x) << "\n";
	std::cout << asin(x) << "\n";
	std::cout << acos(x) << "\n";
	std::cout << atan(x) << "\n";
	std::cout << sinh(x) << "\n";
	std::cout << cosh(x) << "\n";
	std::cout << tanh(x) << "\n";
	std::cout << asinh(x) << "\n";
	std::cout << acosh(x) << "\n";
	std::cout << atanh(x) << "\n";
}
