// test program for interval multiplication

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <limits>

typedef kv::interval<double> itv;

int main()
{
	itv x, y, z;
	itv::base_type inf = std::numeric_limits<itv::base_type>::infinity();

	x = itv(-2, -1);
	y = itv(-2, -1);
	z = x * y;
	if (z.lower() != 1 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, -1);
	y = itv(-2, 3);
	z = x * y;
	if (z.lower() != -6 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, -1);
	y = itv(1, 2);
	z = x * y;
	if (z.lower() != -4 || z.upper() != -1) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, 1);
	y = itv(-2, -1);
	z = x * y;
	if (z.lower() != -2 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, 1);
	y = itv(-2, 2);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, 2);
	y = itv(-1, 2);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, 2);
	y = itv(-2, 1);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-1, 2);
	y = itv(-2, 2);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-2, 1);
	y = itv(1, 2);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 2) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(1, 2);
	y = itv(-2, -1);
	z = x * y;
	if (z.lower() != -4 || z.upper() != -1) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(1, 2);
	y = itv(-2, 1);
	z = x * y;
	if (z.lower() != -4 || z.upper() != 2) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(1, 2);
	y = itv(1, 2);
	z = x * y;
	if (z.lower() != 1 || z.upper() != 4) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, 0);
	y = itv(-inf, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, 0);
	y = itv(-inf, inf);
	z = x * y;
	if (z.lower() != -inf || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, 0);
	y = itv(0, inf);
	z = x * y;
	if (z.lower() != -inf || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, 0);
	y = itv(0, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, inf);
	y = itv(-inf, 0);
	z = x * y;
	if (z.lower() != -inf || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, inf);
	y = itv(-inf, inf);
	z = x * y;
	if (z.lower() != -inf || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, inf);
	y = itv(0, inf);
	z = x * y;
	if (z.lower() != -inf || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(-inf, inf);
	y = itv(0, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, inf);
	y = itv(-inf, 0);
	z = x * y;
	if (z.lower() != -inf || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, inf);
	y = itv(-inf, inf);
	z = x * y;
	if (z.lower() != -inf || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, inf);
	y = itv(0, inf);
	z = x * y;
	if (z.lower() != 0 || z.upper() != inf) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, inf);
	y = itv(0, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, 0);
	y = itv(-inf, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, 0);
	y = itv(-inf, inf);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, 0);
	y = itv(0, inf);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}

	x = itv(0, 0);
	y = itv(0, 0);
	z = x * y;
	if (z.lower() != 0 || z.upper() != 0) {
		std::cout << x << "*" << y << "!=" << z << "\n";
	}
}
