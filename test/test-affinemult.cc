/*
 * test for multiplication of affine arithmetic
 *  use
 *   -DAFFINE_MULT=0 (Stolfi)
 *   -DAFFINE_MULT=1 (better)
 *   -DAFFINE_MULT=2 (best)
 *  and compare the last term of result
 */

#include <kv/affine.hpp>

typedef kv::affine<double> afd;

int main()
{
	afd x, y;

#if AFFINE_SIMPLE >= 1
	x.er = 0;
	y.er = 0;
#endif

/* 
 * below examples are take from:
 *  Shinya MIYAJIMA, Takatomi MIYATA, and Masahide KASHIWAGI:
 *  On the Best Multiplication of the Affine Arithmetic,
 *  IEICE Trans. Vol. J86-A, No.2, pp. 150-159, February, 2003 [In Japanese]
 */ 

	afd::maxnum() = 2;

	x.a.resize(3);
	x.a(0) = -5;
	x.a(1) = 3;
	x.a(2) = -1;

	y.a.resize(3);
	y.a(0) = -1.5;
	y.a(1) = -0.3;
	y.a(2) = -0.2;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";

	afd::maxnum() = 2;

	x.a.resize(3);
	x.a(0) = 3;
	x.a(1) = 1;
	x.a(2) = 0;

	y.a.resize(3);
	y.a(0) = 4;
	y.a(1) = 0;
	y.a(2) = 2;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";

	afd::maxnum() = 3;

	x.a.resize(4);
	x.a(0) = 5;
	x.a(1) = 2;
	x.a(2) = -0.6;
	x.a(3) = 0.4;

	y.a.resize(4);
	y.a(0) = 5;
	y.a(1) = -0.7;
	y.a(2) = -1.1;
	y.a(3) = 0.2;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";

	afd::maxnum() = 4;

	x.a.resize(5);
	x.a(0) = -28;
	x.a(1) = 10;
	x.a(2) = -3;
	x.a(3) = -4;
	x.a(4) = 5;

	y.a.resize(5);
	y.a(0) = -4.5;
	y.a(1) = -1;
	y.a(2) = 0.3;
	y.a(3) = -0.2;
	y.a(4) = 0;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";

	afd::maxnum() = 5;

	x.a.resize(6);
	x.a(0) = -19.5;
	x.a(1) = -4;
	x.a(2) = 0.5;
	x.a(3) = -2;
	x.a(4) = 1;
	x.a(5) = 3;

	y.a.resize(6);
	y.a(0) = -16;
	y.a(1) = -0.3;
	y.a(2) = 1.6;
	y.a(3) = 0.6;
	y.a(4) = 0.1;
	y.a(5) = 1.4;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";

	afd::maxnum() = 6;

	x.a.resize(7);
	x.a(0) = 64;
	x.a(1) = 7;
	x.a(2) = 8;
	x.a(3) = -15;
	x.a(4) = 13;
	x.a(5) = 2;
	x.a(6) = 12;

	y.a.resize(7);
	y.a(0) = 47.5;
	y.a(1) = 20;
	y.a(2) = 1;
	y.a(3) = -14;
	y.a(4) = 3;
	y.a(5) = 6;
	y.a(6) = 1.5;

	std::cout << x << "\n";
	std::cout << y << "\n";
	std::cout << x * y << "\n";
}
