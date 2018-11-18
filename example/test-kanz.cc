#include <boost/timer.hpp>
#include <kv/allsol.hpp>

/*
 * Examples taken from:
 * 
 * Yuchi Kanzawa, Masahide Kashiwagi, Shin'ichi Oishi, Haruyuki Nakamura:
 * An Algorithm of Finding All Solutions with Guaranteed Accuracy
 * for Nonlinear Equations within Finite Steps
 * IEICE transactions A, Vol. J80-A, No.7, pp.1130-1137 (1997.7)
 * (in Japanese)
 * 
 * https://search.ieice.org/bin/summary.php?id=j80-a_7_1130
 */

namespace ub = boost::numeric::ublas;

struct Kanz1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(3);

		static const T a = T(2.);
		static const T b = 0.4 * kv::constants<T>::pi();
		static const T c = T(2.5);
		static const T d = T(1.);
		static const T e = T(0.75);
		static const T f = 2. * kv::constants<T>::pi();
		static const T g = T(1.);
		static const T h = kv::constants<T>::str("0.8");
		static const T k = 2. * kv::constants<T>::pi();

		y(0) = 2. * sin(b * x(0)) * sin(b * x(2)) - x(1);
		y(1) = c - d * x(2) + e * x(1) * sin(f * x(2)) - x(0);
		y(2) = g + h * x(1) * sin(k * x(0)) - x(2);

		return y;
	}
};

struct Kanz2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(4);
		static const T ep = T(1.);
		static const T sigma = T(pow(2., -5));
		static const T omega = T(4.);
		static const T Omega = T(omega * omega);

		y(0) = (9./Omega - 1.) * x(0) - 3. * sigma / omega * x(1) + 9 * ep / Omega * (0.75 * pow(x(0), 3) - 0.75 * pow(x(0), 2) * x(2) + 0.75 * pow(x(1), 2) * x(2) + 0.75 * x(0) * pow(x(1), 2) + 1.5 * x(0) * pow(x(2), 2) + 1.5 * x(0) * pow(x(3), 2) - 1.5 * x(0) * x(1) * x(3));
		y(1) = 3 * sigma / omega * x(0) + (9. / Omega - 1.) * x(1) + 9. * ep / Omega * (0.75 * pow(x(1), 3) + 0.75 * pow(x(0), 2) * x(1) - 0.75 * pow(x(0), 2) * x(3) + 0.75 * pow(x(1), 2) * x(3) + 1.5 * x(1) * pow(x(2), 2) + 1.5 * x(1) * pow(x(3), 2) + 1.5 * x(0) * x(1) * x(2));
		y(2) = (9. / Omega - 9.) * x(2) - 9. * sigma / omega * x(3) + 9. * ep / Omega * (-0.25 * pow(x(0), 3) + 0.75 * pow(x(2), 3) + 1.5 * pow(x(0), 2) * x(2) + 1.5 * pow(x(1), 2) * x(2) + 0.75 * x(0) * pow(x(1), 2) + 0.75 * x(2) * pow(x(3), 2));

		y(3) = (9. / Omega - 9.) * x(3) + 9. * sigma / omega * x(2) - 9. / Omega + 9. * ep / Omega * (0.25 * pow(x(1), 3) + 0.75 * pow(x(3), 3) - 0.75 * pow(x(0), 2) * x(1) + 1.5 * pow(x(0), 2) * x(3) + 1.5 * pow(x(1), 2) * x(3) + 0.75 * pow(x(2), 2) * x(3));

		return y;
	}
};

struct Kanz3 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);
		const T c1 = T("0.9999");
		const T c2 = T("1.0001");
		const T c3 = T("0.0001");

		y(0) = pow(x(0), 3) * x(1) - 0.25 * pow(x(0), 2) + (pow(x(1), 3) - c1 * x(1)) * x(0) - 0.25 * (pow(x(1), 2) - 1.);
		y(1) = pow(x(0), 3) - pow(x(0), 2) * x(1) + (pow(x(1), 2) - c2) * x(0) - (pow(x(1), 3) - x(1) + c3);
		
		return y;
	}
};

struct Kanz4 {
	template <class T> T g1(T x) {
		const T c1 = T("10.5");
		const T c2 = T("11.8");

		return 2.5 * pow(x, 3) - c1 * pow(x, 2) + c2 * x;
	}
	template <class T> T g2(T x) {
		const T c1 = T("0.43");
		const T c2 = T("2.69");
		const T c3 = T("4.56");

		return c1 * pow(x, 3) - c2 * pow(x, 2) + c3 * x;
	}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);
		const T c1 = T("13.3");

		y(0) = 30. - c1 * g1(x(0)) - (x(0) + x(1));
		y(1) = g1(x(0)) - g2(x(1));

		return y;
	}
};

struct Kanz5 {
	int n;

	Kanz5 (int n = 1) : n(n){}
	template <class T> T g(T x) {
		const T c1 = T("10.5");
		const T c2 = T("11.8");

		return 2.5 * pow(x, 3) - c1 * pow(x, 2) + c2 * x;
	}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int i;
		ub::vector<T> y(n);
		T s(0.);

		for (i=0; i<n; i++) {
			s += x(i);
		}

		for (i=0; i<n; i++) {
			y(i) = g(x(i)) + s - T(i+1);
		}

		return y;
	}
};

struct Kanz6 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(3);
		const T c = T("3.84");

		y(0) = x(1) - c * x(0) * (1. - x(0));
		y(1) = x(2) - c * x(1) * (1. - x(1));
		y(2) = x(0) - c * x(2) * (1. - x(2));

		return y;
	}
};


int main()
{
	int i, n;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	std::cout.precision(17);


	I.resize(3);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-5., 5.);

	t.restart();
	allsol(Kanz1(), I);
	std::cout << t.elapsed() << " sec\n";

	I.resize(4);
	for (i=0; i<2; i++) I(i) = kv::interval<double>(-3., 3.);
	for (i=2; i<4; i++) I(i) = kv::interval<double>(-0.3, 0.3);

	t.restart();
	allsol(Kanz2(), I);
	std::cout << t.elapsed() << " sec\n";

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-1., 1.);

	t.restart();
	allsol(Kanz3(), I);
	std::cout << t.elapsed() << " sec\n";

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-5., 5.);

	t.restart();
	allsol(Kanz4(), I);
	std::cout << t.elapsed() << " sec\n";

	for (n=3; n<=6; n++) {
		I.resize(n);
		for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-2., 2.);

		t.restart();
		allsol(Kanz5(n), I);
		std::cout << t.elapsed() << " sec\n";
	}

	I.resize(3);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(0., 1.);

	t.restart();
	allsol(Kanz6(), I);
	std::cout << t.elapsed() << " sec\n";
}
