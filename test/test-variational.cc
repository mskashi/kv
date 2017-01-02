#include <kv/ode.hpp>
#include <kv/ode-autodif.hpp>
#include <kv/autodif.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

namespace kv {

template <class F> struct VariationalEq {
	F f;

	VariationalEq(F func) : f(func) {}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
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
	ub::vector<itv> x;
	ub::vector< kv::autodif<itv> > dx;
	int r;
	itv end;

	std::cout.precision(17);

	// standard calculation

	Lorenz f;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;

	r = kv::ode(f, x, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// using automatic differentiation

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	dx = kv::autodif<itv>::init(x);
	end = 1.;

	r = kv::ode(f, dx, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else {
		std::cout << dx << "\n";
		std::cout << end << "\n";
	}

	// using "variational equation maker"

	kv::VariationalEq<Lorenz> g(f);

	x.resize(3 + 3 * 3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	x(3) = 1.; x(4) = 0.; x(5) = 0.;
	x(6) = 0.; x(7) = 1.; x(8) = 0.;
	x(9) = 0.; x(10) = 0.; x(11) = 1.;
	end = 1.;

	r = kv::ode(g, x, itv(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
