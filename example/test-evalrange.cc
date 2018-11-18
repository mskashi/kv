#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/affine.hpp>
#include <kv/autodif.hpp>
#include <kv/interval-vector.hpp>

/*
 * comparison of interval, affine, mean value form
 * in evaluation of a 5 variable function.
 * The function is taken from http://ci.nii.ac.jp/naid/110003291707
 */

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;

struct Func {
	template <class T> T operator()(const ub::vector<T>& t) {
		T f1, f2, y;

		f1 = (-10. - 13.2 * t[0] + 0.5 * t[1] - 0.4 * t[2] + 7.3 * t[3] + 4.4 * t[4]);
		f2 = (0.8 + 10.5 * t[0] + 5.9 * t[1] + 7.5 * t[2] - 14.5 * t[3] + 1.3 * t[4]);

		y = f1 * f2 - ((-113.7) * f2 + (-49.8) * f1) + 5660.;
		return y;
	}
};

int main()
{
	int i;
	ub::vector<itv> I;
	Func f;

	std::cout.precision(17);

	I.resize(5);

	I[0] = itv(8.7, 8.8);
	I[1] = itv(-9.4, -9.3);
	I[2] = itv(-4.6, -4.5);
	I[3] = itv(3.5, 3.6);
	I[4] = itv(-2.9, -2.8);


	// interval arithmetic

	std::cout << "interval: " << f(I) << "\n";

	// affine arithmetic

	// ub::vector< kv::affine<double> > a;
	ub::vector< kv::affine<itv::base_type> > a;

	a.resize(5);
	for (i=0; i<5; i++) a(i) = I(i);
	std::cout << "affine: " << to_interval(f(a)) << "\n";

	// mean value form

	ub::vector<itv> c, fdi;
	itv fc, fi;

	c = mid(I);
	fc = f(c);
	kv::autodif<itv>::split(f(kv::autodif<itv>::init(I)), fi, fdi);

	std::cout << "mvf: " << fc + inner_prod(fdi, I - c) << "\n";
}
