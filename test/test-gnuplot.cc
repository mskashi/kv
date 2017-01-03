#include <kv/gnuplot.hpp>

int main()
{
	kv::gnuplot g;

	// initialize
	g.open();

	// initialize with the path of gnuplot exectable
	// g.open("/usr/bin/gnuplot");

	// example for generating picture file
	// g.open("/usr/bin/gnuplot", "postscript eps color", "aaa.eps");

	// example for saving gnuplot command sequence using tee
	// g.open("tee aaa.dat | /usr/bin/gnuplot");

	// set drawing range
	g.screen(0, 0, 10, 10);

	g.line(1,1,3,4);
	g.point(4,4);
	g.ellipsef(8,2,2,1);
	g.rect(5,5,6,7);
	g.rectf(8,8,9,9);

	getchar();

	// finish drawing
	g.close();
}
