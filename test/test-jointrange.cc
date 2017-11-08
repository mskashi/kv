#include <kv/affine.hpp>
#include <kv/jointrange.hpp>

typedef kv::affine<double> afd;

int main()
{
	kv::matplotlib g;
	g.open();
	g.screen(-3,-3,3,3);
	g.send_command("set size ratio -1");

	afd x, y;

	afd::maxnum() = 3;

	x.a.resize(4);
	x.a(0) = 0;
	x.a(1) = 1;
	x.a(2) = -1;
	x.a(3) = 0.1;

	y.a.resize(4);
	y.a(0) = 0;
	y.a(1) = 1;
	y.a(2) = 1;
	y.a(3) = 0.2;

#if AFFINE_SIMPLE >= 1
	x.er = 0.1;
	y.er = 0.1;
#endif

	std::cout << x << "\n";
	std::cout << y << "\n";

	kv::jointrange(x, y, g);

	getchar();

	g.close();
}
