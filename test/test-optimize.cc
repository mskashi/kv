#include <kv/optimize.hpp>

namespace ub = boost::numeric::ublas;

struct Func {
	template <class T> T operator() (const ub::vector<T>& x){
		return 1/(pow(x(0)-2,2) + pow(x(1)-3,2)+1) - 1/(pow(x(0)-5,2)+pow(x(1)-6,2)+1);
	}
};

struct Func2 {
	template <class T> T operator() (const T& x){
		return sin(x) + log(x);
	}
};

// example found in
// Numerical Optimizer/Global by NTT Data
// http://www.msi.co.jp/nuopt/products/derivation/global/index.html
struct NTTData {
	template <class T> T operator() (const ub::vector<T>& in){
		T x, y, z, w;

		x = in(0);
		y = in(1);
		z = in(2);
		w = in(3);

		return 100. * (y - x*x) * (y - x*x) + (1. - x ) * (1. - x)
		    + 90. * (z - w * w) * (z - w * w) + (1. - w) * (1. - w)
		    + 10.1 * (y - 1.) * (y - 1.) + 10.1 * (z - 1.) * (z - 1.)
		    + 19.8 * (y - 1.) * (z - 1.);
	}
};

struct Rosenbrock {
	template <class T> T operator() (const ub::vector<T>& x){
		T tmp, tmp2;
		tmp = 1. - x(0);
		tmp2 = x(1) - x(0) * x(0);
		return tmp * tmp + 100. * tmp2 * tmp2;
	}
};

// taken from: http://ci.nii.ac.jp/naid/110002373498

struct Levy {
	template <class T> T operator() (const ub::vector<T>& x){
		T tmp, tmp2;
		int i;

		tmp = 0;
		for (i=1; i<=5; i++) {
			tmp += i * cos((i-1)*x(0) + i);
		}
		tmp2 = 0;
		for (i=1; i<=5; i++) {
			tmp2 += i * cos((i+1)*x(1) + i);
		}
		
		return tmp * tmp2;
	}
};

/*
 * Griewank Function
    taken from:
 *  https://www.sfu.ca/~ssurjano/griewank.html
 *  http://www.autodiff.org/ad16/Invited/Rump_Reliable.pdf
 */
struct Griewank {
	int n;
	Griewank(int n): n(n) {
	}
	template <class T> T operator() (const ub::vector<T>& x){
		int i;
		T s1, s2;
		s1 = 0;
		for (i=0; i<n; i++) {
			s1 += x(i) * x(i);
		}
		s2 = 1;
		for (i=0; i<n; i++) {
			s2 *= cos(x(i) / sqrt(T(i + 1)));
		}
		return 1 + s1 / 4000. - s2;
	}
};


typedef kv::interval<double> itv;

int main()
{
	int i, n;
	ub::vector<itv> I;
	std::list< ub::vector<itv> > result;
	std::list< ub::vector<itv> >::iterator p;
	std::list<itv> result2;
	std::list<itv>::iterator p2;
	itv v;

	std::cout.precision(17);

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = itv(0., 10.);

	std::cout << "calculate intervals minimizing Func with partitioned intervals smaller than 1e-5\n";
	result = minimize(I, Func(), 1e-5);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}

	std::cout << "calculate intervals minimizing Func with 10000 partitioned intervals\n";
	result = minimize(I, Func(), 10000);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}

	std::cout << "calculate minimum value of Func\n";
	v = minimize_value(I, Func(), 10000);
	std::cout << v << "\n";

	std::cout << "calculate intervals maximizing Func\n";
	result = maximize(I, Func(), 10000);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}

	std::cout << "calculate maximum value of Func\n";
	v = maximize_value(I, Func(), 10000);
	std::cout << v << "\n";

	std::cout << "calculate intervals minimizing Func2\n";
	result2 = minimize(itv(1,10), Func2(), 10000);
	p2 = result2.begin();
	while (p2 != result2.end()) {
		std::cout << *(p2++) << "\n";
	}

	std::cout << "calculate minimum value of Func2\n";
	v = minimize_value(itv(1,10), Func2(), 10000);
	std::cout << v << "\n";

	std::cout << "calculate intervals maximizing Func2\n";
	result2 = maximize(itv(1,10), Func2(), 10000);
	p2 = result2.begin();
	while (p2 != result2.end()) {
		std::cout << *(p2++) << "\n";
	}

	std::cout << "calculate maximum value of Func2\n";
	v = maximize_value(itv(1,10), Func2(), 10000);
	std::cout << v << "\n";


	std::cout << "NTTData\n";
	I.resize(4);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);
	result = minimize(I, NTTData(), 1e-3);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}

	std::cout << "Levy\n";
	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = itv(0., 10.);
	result = minimize(I, Levy(), 1e-5);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}

	std::cout << "Griewank\n";
	n = 5;
	I.resize(n);
	for (i=0; i<I.size(); i++) I(i) = itv(-600., 600.);
	result = minimize(I, Griewank(n), 1e-5);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}
}
