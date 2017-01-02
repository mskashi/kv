#include "interval.hpp"
#include "rdouble.hpp"

#include "psa.hpp"
#include "psa-plot.hpp"

typedef kv::interval<double> itvd;


int main()
{
	kv::psa<itvd> c;
	kv::gnuplot g;
	g.open("/usr/bin/gnuplot", "x11");
	g.screen(-0.1,-5,0.3,5);

	c.v.resize(3);
	c.v(0) = 1.;
	c.v(1) = 3.;
	c.v(2) = 2.;

	kv::psa<itvd>::domain() = itvd(0., 0.2);

	std::cout << c << "\n";
	kv::psa_plot(c, 0., g);
	getchar();

	std::cout << c*c << "\n";
	kv::psa_plot(c*c, 0., g);

	getchar();

	c = c/c;
	std::cout << c << "\n";
	kv::psa_plot(c, 0., g);

	getchar();
}
