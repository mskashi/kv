#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <kv/psa.hpp>
#include <kv/psa-plot.hpp>

typedef kv::interval<double> itv;


int main()
{
	kv::psa<itv> c;
	kv::matplotlib g;
	g.open();
	g.screen(-0.1,-5,0.3,5);

	c.v.resize(3);
	c.v(0) = 1.;
	c.v(1) = 3.;
	c.v(2) = 2.;

	kv::psa<itv>::domain() = itv(0., 0.2);

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
