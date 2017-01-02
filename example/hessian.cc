/*
 * Calculating Hessian of f:R^n -> R using nested autodif
 */

#include <kv/autodif.hpp>

namespace ub = boost::numeric::ublas;

int main()
{
	ub::vector< kv::autodif< kv::autodif<double> > > ax;
	kv::autodif< kv::autodif<double> > ay;
	kv::autodif<double> tmpv;
	ub::vector< kv::autodif<double> > tmpd;
	ub::vector<double> tmpdv;
	ub::matrix<double> tmpdd;

	ub::vector<double> x;
	int i, j;

	x.resize(2);
	x(0) = 2.;
	x(1) = 3.;

	ax = kv::autodif< kv::autodif<double> >::init(kv::autodif<double>::init(x));

	ay = 1 / (pow(ax(0), 3) + pow(ax(1), 2));

	std::cout << ay << "\n";

	kv::autodif< kv::autodif<double> >::split(ay, tmpv, tmpd);
	kv::autodif<double>::split(tmpd, tmpdv, tmpdd);

	std::cout << "f:\n";
	std::cout << tmpv.v << "\n";
	std::cout << "f':\n";
	std::cout << tmpv.d << "\n";
	std::cout << "f':\n";
	std::cout << tmpdv << "\n";
	std::cout << "f'':\n";
	std::cout << tmpdd << "\n";

	/*
		(check with mathematica)
		x=2; y=3; 1/(x^3+y^2)
		x=2; y=3; D[1/(x^3+y^2),x]
		x=2; y=3; D[1/(x^3+y^2),y]
		x=2; y=3; D[D[1/(x^3+y^2),x],x]
		x=2; y=3; D[D[1/(x^3+y^2),x],y]
		x=2; y=3; D[D[1/(x^3+y^2),y],x]
		x=2; y=3; D[D[1/(x^3+y^2),y],y]
	*/
}
