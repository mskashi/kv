#include "ode.hpp"
#include "autodif.hpp"

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};

class Lorenz {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

namespace kv {

template <class F> class VariationalEq {
	public:
	F f;

	VariationalEq(F func) : f(func) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		int s = x.size();
		int s2 = (int)std::floor((-1.+std::sqrt(1. + 4.*s))/2.+0.5);

		ub::vector<T> x1(s2);
		ub::matrix<T> x2(s2, s2);
		ub::vector<T> y(s);

		ub::vector< autodif<T> > r;
		ub::vector<T> rv;
		ub::matrix<T> rm;

		int i, j, k;

		k = 0;
		for (i=0; i<s2; i++) {
			x1(i) = x(k++);
		}
		for (i=0; i<s2; i++) {
			for (j=0; j<s2; j++) {
				x2(i, j) = x(k++);
			}
		}

		autodif<T>::split(f(autodif<T>::init(x1), autodif<T>(t)), rv, rm);

		rm = prod(rm, x2);

		k = 0;
		for (i=0; i<s2; i++) {
			y(k++) = rv(i);
		}
		for (i=0; i<s2; i++) {
			for (j=0; j<s2; j++) {
				y(k++) = rm(i, j);
			}
		}

		return y;
	}
};

} // namespace kv



int main()
{
	ub::vector<itvd> x;
	bool r;

	itvd end;

	end = std::numeric_limits<double>::infinity();

	// x.resize(2);
	// x(0) = 1.; x(1) = 0.;

	// x.resize(3);
	// x(0) = 15.; x(1) = 15.; x(2) = 36.;
	x.resize(3 + 3 * 3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	x(3) = 1.; x(4) = 0.; x(5) = 0.;
	x(6) = 0.; x(7) = 1.; x(8) = 0.;
	x(9) = 0.; x(10) = 0.; x(11) = 1.;

	std::cout.precision(17);

	Lorenz f;
	kv::VariationalEq<Lorenz> g(f);

	// r = kv::ode(Func(), x, itvd(0.), itvd(1.0), 5);
	r = kv::ode(g, x, itvd(0.), end, 12);

	if (!r) std::cout << "No Solution\n";
	else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
