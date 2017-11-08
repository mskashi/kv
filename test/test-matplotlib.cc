#include <kv/matplotlib.hpp>

int main()
{
	kv::matplotlib g;

	// initialize
	g.open();

	// set drawing range
	g.screen(0, 0, 10, 10);
	// aspect ratio = 1
	// g.screen(0, 0, 10, 10, true);

	g.line(1,1,3,4);
	g.line(1,2,3,5, "red");

	g.rect(4,6,5,8);
	g.rect(6,6,7,8, "green");
	g.rect(8,6,9,8, "black", "red");

	g.point(4,2);
	g.ellipse(8,2,2,1);
	g.circle(2,8,2);

	double xs[5] = {6., 7., 6., 5., 5.};
	double ys[5] = {4., 5., 6., 6., 5.};
	g.polygon(xs, ys, 5, "black", "yellow");

	g.line(1,3,3,6, "green", "alpha=0.2");
	g.line(1,4,4,5, "black", "alpha=0.5, linestyle='--'");
	g.point(4, 3, "red", "s=100");

	g.save("test.pdf");

	getchar();

	// finish drawing
	g.close();
}
