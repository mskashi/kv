#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/affine.hpp>
#include <kv/autodif.hpp>
#include <kv/interval-vector.hpp>

/*
 * http://ci.nii.ac.jp/naid/110003291707 ¤Ë¤¢¤Ã¤¿5ÊÑ¿ô´Ø¿ô¤ÎÁü¤ÎÉ¾²Á¤ò
 * ¶è´Ö±é»»¡¢affine, Ê¿¶ÑÃÍ·Á¼°¤ÇÈæ³Ó
 */

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itvd;

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
	ub::vector<itvd> I;
	Func f;

	std::cout.precision(17);

	I.resize(5);

	I[0] = itvd(8.7, 8.8);
	I[1] = itvd(-9.4, -9.3);
	I[2] = itvd(-4.6, -4.5);
	I[3] = itvd(3.5, 3.6);
	I[4] = itvd(-2.9, -2.8);

	// interval arithmetic
	std::cout << "interval: " << f(I) << "\n";


	ub::vector< kv::affine<double> > a;

	a.resize(5);

	for (i=0; i<5; i++) a(i) = I(i);

	// affine arithmetic
	std::cout << "affine: " << to_interval(f(a)) << "\n";


	ub::vector<itvd> c, fdi;
	itvd fc, fi;

	c = mid(I);
	fc = f(c);
	kv::autodif<itvd>::split(f(kv::autodif<itvd>::init(I)), fi, fdi);

	// mean value form
	std::cout << "mvf: " << fc + inner_prod(fdi, I - c) << "\n";
}
