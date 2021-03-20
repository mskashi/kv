#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/autodif.hpp>

namespace ub = boost::numeric::ublas;

// f: R->R
kv::autodif<double> testfunc1(const kv::autodif<double>& x) {
	kv::autodif<double> tmp;

	tmp = x * x - 1.;
	return cos(tmp + sin(tmp));
}

// f: R->R (use template)
template <class T> T testfunc2(const T& x) {
	T tmp;

	tmp = x * x - 1.;
	return cos(tmp + sin(tmp));
}

// f: R^n->R
template <class T> T testfunc3(const ub::vector<T>& x) {
	return atan(x(0) * x(0) + 2. * x(1) * x(1) - 1.);
}

// f: R->R^n
template <class T> ub::vector<T> testfunc4(const T& x) {
	ub::vector<T> y(2);

	y(0) = sin(x);
	y(1) = cos(x);

	return y;
}

// f: R^n->R^n
template <class T> ub::vector<T> testfunc5(const ub::vector<T>& x) {
	ub::vector<T> y(2);

	y(0) = 2. * x(0) * x(0) * x(1) - 1.;
	y(1) = x(0) + 0.5 * x(1) * x(1) - 2.;

	return y;
}


int main()
{
	double d1, d2;
	kv::autodif<double> a1, a2, a3;
	ub::vector<double> v1, v2;
	ub::vector< kv::autodif<double> > va1, va2, va3, va4, va5;
	ub::matrix<double> m;
	int i;

	//
	// f: R->R (testfunc1)
	//

	// helper function to initialize autodif variable with single variable
	// same as " a1.v = 1.5; a1.d(0) = 1.; "
	a1 = kv::autodif<double>::init(1.5);

	a2 = testfunc1(a1);

	// helper function to separate autodif variable into value and gradient
	// same as " d1 = a2.v; d2 = a2.d(0); "
	kv::autodif<double>::split(a2, d1, d2);

	std::cout << "testfunc1\n";
	std::cout << d1 << "\n"; // f(1.5)
	std::cout << d2 << "\n"; // f'(1.5)

	//
	// f: R->R (testfunc2: use template)
	//

	a1 = kv::autodif<double>::init(1.5);
	a2 = testfunc2(a1);

	kv::autodif<double>::split(a2, d1, d2);

	std::cout << "testfunc2\n";
	std::cout << d1 << "\n"; // f(1.5)
	std::cout << d2 << "\n"; // f'(1.5)

	d1 = 1.5;
	d2 = testfunc2(d1); // double can be also passed
	std::cout << d2 << "\n"; // f(1.5)

	//
	// f: R^n->R
	//

	v1.resize(2);
	v1(0) = 5.; v1(1) = 6.;

	// helper function to initialize autodif vector with vector variable
	// same as " va1(0).v = 5.; va1(0).d(0) = 1.; va1(0).d(1) = 0,; "
	//         " va1(1).v = 6.; va1(1).d(0) = 0.; va1(1).d(1) = 1,; "
	va1 = kv::autodif<double>::init(v1);

	a1 = testfunc3(va1);

	// helper function to separate autodif variable into value and gradient
	// same as " y = a1.v; v2(0) = a1.d(0); v2(1) = a1.d(1) "
	kv::autodif<double>::split(a1, d1, v2);

	std::cout << "testfunc3\n";
	std::cout << d1 << "\n"; // f(5, 6)
	std::cout << v2 << "\n"; // gradient

	//
	// f: R->R^n
	//

	a1 = kv::autodif<double>::init(1.5);

	va1 = testfunc4(a1);

	// helper function to separate autodif variable into value and gradient
	// same as " v1(0) = va1(0).v; v1(1) = va1(1).v; "
	//         " v2(0) = va1(0).d(0); v2(1) = va1(1).d(0); "
	kv::autodif<double>::split(va1, v1, v2);

	std::cout << "testfunc4\n";
	std::cout << v1 << "\n"; // f(1.5)
	std::cout << v2 << "\n"; // f'(1.5)

	//
	// f: R^n->R^m
	//

	v1(0) = 5.; v1(1) = 6.;

	va1 = kv::autodif<double>::init(v1);

	va2 = testfunc5(va1);

	// helper function to separate autodif vector into vector value and jacobian
	// same as " v2(0) = va2(0).v; m(0)(0) = va2(0).d(0); m(0)(1) = va2(0).d(1); "
	//         " v2(1) = va2(1).v; m(1)(0) = va2(1).d(0); m(1)(1) = va2(1).d(1); "
	kv::autodif<double>::split(va2, v2, m);

	std::cout << "testfunc5\n";
	std::cout << v2 << "\n"; // f(5, 6)
	std::cout << m << "\n"; // Jacobian matrix

	//
	// operations for autodif
	//

	v1(0) = 0.7; v1(1) = 0.8;
	va1 = kv::autodif<double>::init(v1);
	a1 = va1(0);
	a2 = va1(1);

	std::cout << "autodif operations\n";

	std::cout << a1 << "\n";
	std::cout << a2 << "\n";

	a3 = a1 + a2; std::cout << a3 << "\n";
	a3 = a1 - a2; std::cout << a3 << "\n";
	a3 = a1 * a2; std::cout << a3 << "\n";
	a3 = a1 / a2; std::cout << a3 << "\n";
	a3 = pow(a1, 3); std::cout << a3 << "\n";
	a3 = pow(a1, a2); std::cout << a3 << "\n";
	a2 = sqrt(a1); std::cout << a2 << "\n";
	a2 = exp(a1); std::cout << a2 << "\n";
	a2 = log(a1); std::cout << a2 << "\n";
	a2 = sin(a1); std::cout << a2 << "\n";
	a2 = cos(a1); std::cout << a2 << "\n";
	a2 = tan(a1); std::cout << a2 << "\n";
	a2 = sinh(a1); std::cout << a2 << "\n";
	a2 = cosh(a1); std::cout << a2 << "\n";
	a2 = tanh(a1); std::cout << a2 << "\n";
	a2 = asin(a1); std::cout << a2 << "\n";
	a2 = acos(a1); std::cout << a2 << "\n";
	a2 = atan(a1); std::cout << a2 << "\n";
	// asinh, acosh, atanh may not compile for autodif<double>
	// because c++03 does not have asinh, acosh, atanh
	// a2 = asinh(a1); std::cout << a2 << "\n";
	// a2 = acosh(a1+1.); std::cout << a2 << "\n";
	// a2 = atanh(a1); std::cout << a2 << "\n";

	//
	// compress / expand
	//

	v1.resize(20);
	for (i=0; i<20; i++) v1(i) = i + 2.;
	va1 = kv::autodif<double>::init(v1);
	// use only a few of the many variables
	va2.resize(2);
	va2(0) = va1(10);
	va2(1) = va1(11);

	// usual way
	std::cout << "usual input\n";
	std::cout << va2 << "\n";
	va3 = testfunc5(va2);
	std::cout << "usual output\n";
	std::cout << va3 << "\n";

	// using compress/expand to save calculation time
	ub::matrix<double> save;
	va3 = kv::autodif<double>::compress(va2, save);
	std::cout << "compressed input\n";
	std::cout << va3 << "\n";
	va4 = testfunc5(va3);
	std::cout << "compressed output\n";
	std::cout << va4 << "\n";
	va5 = kv::autodif<double>::expand(va4, save);
	std::cout << "expanded output\n";
	std::cout << va5 << "\n";
}
