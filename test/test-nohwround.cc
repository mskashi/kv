/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

/*
 * test program for rounded math emulation
 *  compare add, sub, mul, div and sqrt under
 *  - changing rounding mode
 *  - emulating rounded math using twosum and twoproduct
 */

#include <fenv.h>
#include <boost/random.hpp>
#include <ctime>

#ifndef NT
#define NT 10000000
#endif

struct hwround {
	public:

	static void roundnear() {
		fesetround(FE_TONEAREST);
	}

	static void rounddown() {
		fesetround(FE_DOWNWARD);
	}

	static void roundup() {
		fesetround(FE_UPWARD);
	}

	static void roundchop() {
		fesetround(FE_TOWARDZERO);
	}

	static double add_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 + y1;
		roundnear();
		return r;
	}

	static double add_down(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 + y1;
		roundnear();
		return r;
	}

	static double sub_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 - y1;
		roundnear();
		return r;
	}

	static double sub_down(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 - y1;
		roundnear();
		return r;
	}

	static double mul_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 * y1;
		roundnear();
		return r;
	}

	static double mul_down(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 * y1;
		roundnear();
		return r;
	}

	static double div_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 / y1;
		roundnear();
		return r;
	}

	static double div_down(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 / y1;
		roundnear();
		return r;
	}

	static double sqrt_up(const double& x) {
		volatile double r, x1 = x;
		roundup();
		r = sqrt(x1);
		roundnear();
		return r;
	}

	static double sqrt_down(const double& x) {
		volatile double r, x1 = x;
		rounddown();
		r = sqrt(x1);
		roundnear();
		return r;
	}
};

struct nohwround {

	static void fasttwosum(const double& a, const double& b, double& x, double& y) {
		double tmp;
		x = a + b;
		tmp = x - a;
		y = b - tmp;
	}

	#if 0
	static void twosum(const double& a, const double& b, double& x, double& y) {
		double tmp;
		x = a + b;
		tmp = x - a;
		if (std::fabs(tmp) == std::numeric_limits<double>::infinity()) {
			tmp = x * 0.5 - a * 0.5;
			y = ((a * 0.5 - (x * 0.5 - tmp)) + (b * 0.5 - tmp)) * 2.;
		} else {
			y = (a - (x - tmp)) + (b - tmp);
		}
	}
	#endif

	static void twosum(const double& a, const double& b, double& x, double& y) {
		double tmp;

		x = a + b;
		if (std::fabs(a) > std::fabs(b)) {
			tmp = x - a;
			y = b - tmp;
		} else {
			tmp = x - b;
			y = a - tmp;
		}
	}

	static void split(const double& a, double& x, double& y) {
		static const double sigma = ldexp(1., 27) + 1.;
		double tmp;

		tmp = a * sigma;
		x = tmp - (tmp - a);
		y = a - x;
	}

	static void twoproduct(const double& a, const double& b, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);
		static const double th2 = ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;

		x = a * b;
		#if 0
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
		#endif
		if (std::fabs(a) > th) {
			na = a * c1;
			nb = b * c2;
		} else if (std::fabs(b) > th) {
			na = a * c2;
			nb = b * c1;
		} else {
			na = a;
			nb = b;
		}
		split(na, a1, a2);
		split(nb, b1, b2);
		if (std::fabs(x) > th2) {
			y = a2 * b2 - ((((x * 0.5) - (a1 * 0.5)  * b1) * 2. - a2 * b1) - a1 * b2);
		} else {
			y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
		}
	}

	// succ and pred by Rump

	static double succ(const double& x) {
		static const double th1 = ldexp(1., -969);
		static const double th2 = ldexp(1., -1021);
		static const double c1 = ldexp(1., -53) + ldexp(1., -105);
		static const double c2 = ldexp(1., -1074);
		static const double c3 = ldexp(1., 53);
		static const double c4 = ldexp(1., -53);

		double a, c, e;

		a = std::fabs(x);
		if (a >= th1) return x + a * c1;
		if (a < th2) return x + c2;
		c = c3 * x;
		e = c1 * std::fabs(c);
		return (c + e) * c4;
	}

	static double pred(const double& x) {
		static const double th1 = ldexp(1., -969);
		static const double th2 = ldexp(1., -1021);
		static const double c1 = ldexp(1., -53) + ldexp(1., -105);
		static const double c2 = ldexp(1., -1074);
		static const double c3 = ldexp(1., 53);
		static const double c4 = ldexp(1., -53);

		double a, c, e;

		a = std::fabs(x);
		if (a >= th1) return x - a * c1;
		if (a < th2) return x - c2;
		c = c3 * x;
		e = c1 * std::fabs(c);
		return (c - e) * c4;
	}


	static double add_up(const double& x, const double& y) {
		double r, r2;

		twosum(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return r;
		} else if (r == -std::numeric_limits<double>::infinity()) {
			if (x == -std::numeric_limits<double>::infinity() || y == -std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return -(std::numeric_limits<double>::max)();
			}
		}

		if (r2 > 0.) {
			return succ(r);
		}

		return r;
	}

	static double add_down(const double& x, const double& y) {
		double r, r2;

		twosum(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			if (x == std::numeric_limits<double>::infinity() || y == std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return (std::numeric_limits<double>::max)();
			}
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}

		if (r2 < 0.) {
			return pred(r);
		}

		return r;
	}

	static double sub_up(const double& x, const double& y) {
		double r, r2;

		twosum(x, -y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return r;
		} else if (r == -std::numeric_limits<double>::infinity()) {
			if (x == -std::numeric_limits<double>::infinity() || y == std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return -(std::numeric_limits<double>::max)();
			}
		}

		if (r2 > 0.) {
			return succ(r);
		}

		return r;
	}

	static double sub_down(const double& x, const double& y) {
		double r, r2;

		twosum(x, -y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			if (x == std::numeric_limits<double>::infinity() || y == -std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return (std::numeric_limits<double>::max)();
			}
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}

		if (r2 < 0.) {
			return pred(r);
		}

		return r;
	}

	static double mul_up(const double& x, const double& y) {
		double r, r2;
		double x1, y1;
		double s, s2, t;
		static const double th = ldexp(1., -969); // -1074 + 106 - 1
		static const double c = ldexp(1., 537); // 1074 / 2

		// if (x == 0. || y == 0.) return x * y;

		twoproduct(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return r;
		} else if (r == -std::numeric_limits<double>::infinity()) {
			if (std::fabs(x) == std::numeric_limits<double>::infinity() || std::fabs(y) == std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return -(std::numeric_limits<double>::max)();
			}
		}

		if (fabs(r) >= th) {
			if (r2 > 0.) return succ(r);
			return r;
		} else {
			twoproduct(x * c, y * c, s, s2);
			t = (r * c) * c;
			if ( t < s || (t == s && s2 > 0.)) {
				return succ(r);
			}
			return r;
		}
	}

	static double mul_down(const double& x, const double& y) {
		double r, r2;
		double x1, y1;
		double s, s2, t;
		static const double th = ldexp(1., -969); // -1074 + 106 - 1
		static const double c = ldexp(1., 537); // 1074 / 2

		// if (x == 0. || y == 0.) return x * y;

		twoproduct(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			if (std::fabs(x) == std::numeric_limits<double>::infinity() || std::fabs(y) == std::numeric_limits<double>::infinity()) {
				return r;
			} else {
				return (std::numeric_limits<double>::max)();
			}
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}

		if (fabs(r) >= th) {
			if (r2 < 0.) return pred(r);
			return r;
		} else {
			twoproduct(x * c, y * c, s, s2);
			t = (r * c) * c;
			if ( t > s || (t == s && s2 < 0.)) {
				return pred(r);
			}
			return r;
		}
	}

	static double div_up(const double& x, const double& y) {
		double r, r2;
		double xn, yn, d;
		static const double th1 = ldexp(1., -969); // -1074 + 106 - 1
		static const double th2 = ldexp(1., 918); // 1023 - 105
		static const double c1 = ldexp(1., 105); // -969 - (-1074)
		static const double c2 = ldexp(1., -1074);

		if (x == 0. || y == 0. || std::fabs(x) == std::numeric_limits<double>::infinity() || std::fabs(y) == std::numeric_limits<double>::infinity() || x != x  || y != y) {
			return x / y;
		}

		if (y < 0.) {
			xn = -x;
			yn = -y;
		} else {
			xn = x;
			yn = y;
		}

		if (fabs(xn) < th1) {
			if (fabs(yn) < th2) {
				xn *= c1;
				yn *= c1;
			} else {
				if (xn < 0.) return 0.;
				else return c2;
			}
		}

		d = xn / yn;

		if (d == std::numeric_limits<double>::infinity()) {
			return d;
		} else if (d == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<double>::max)();
		}

		twoproduct(d, yn, r, r2);
		if ( r < xn || ((r == xn) && r2 < 0.)) {
			return succ(d);
		}
		return d;
	}

	static double div_down(const double& x, const double& y) {
		double r, r2;
		double xn, yn, d;
		static const double th1 = ldexp(1., -969); // -1074 + 106 - 1
		static const double th2 = ldexp(1., 918); // 1023 - 105
		static const double c1 = ldexp(1., 105); // -969 - (-1074)
		static const double c2 = ldexp(1., -1074);

		if (x == 0. || y == 0. || std::fabs(x) == std::numeric_limits<double>::infinity() || std::fabs(y) == std::numeric_limits<double>::infinity() || x != x  || y != y) {
			return x / y;
		}

		if (y < 0.) {
			xn = -x;
			yn = -y;
		} else {
			xn = x;
			yn = y;
		}

		if (fabs(xn) < th1) {
			if (fabs(yn) < th2) {
				xn *= c1;
				yn *= c1;
			} else {
				if (xn < 0.) return -c2;
				else return 0.;
			}
		}

		d = xn / yn;

		if (d == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<double>::max)();
		} else if (d == -std::numeric_limits<double>::infinity()) {
			return d;
		}

		twoproduct(d, yn, r, r2);
		if ( r > xn || ((r == xn) && r2 > 0.)) {
			return pred(d);
		}
		return d;
	}

	static double sqrt_up(const double& x) {
		double r, r2, d;
		static const double th1 = ldexp(1., -969); // -1074 + 106 - 1
		static const double c1 = ldexp(1., 106); // -969 - (-1074) + 1
		static const double c2 = ldexp(1., 53); // sqrt(c1)

		d = sqrt(x);

		if (x < th1) {
			double d2, x2;
			x2 = x * c1;
			d2 = d * c2;
			twoproduct(d2, d2, r, r2);
			if ( r < x2 || (r == x2 && r2 < 0.)) {
				return succ(d);
			}
			return d;
		}

		twoproduct(d, d, r, r2);
		if ( r < x || (r == x && r2 < 0.)) {
			return succ(d);
		}
		return d;
	}

	static double sqrt_down(const double& x) {
		double r, r2, d;
		static const double th1 = ldexp(1., -969); // -1074 + 106 - 1
		static const double c1 = ldexp(1., 106); // -969 - (-1074) + 1
		static const double c2 = ldexp(1., 53); // sqrt(c1)

		d = sqrt(x);

		if (x < th1) {
			double d2, x2;
			x2 = x * c1;
			d2 = d * c2;
			twoproduct(d2, d2, r, r2);
			if ( r > x2 || (r == x2 && r2 > 0.)) {
				return pred(d);
			}
			return d;
		}

		twoproduct(d, d, r, r2);
		if ( r > x || (r == x && r2 > 0.)) {
			return pred(d);
		}
		return d;
	}

};

bool samedouble(double x, double y)
{
	// return *((unsigned long long *)(&x)) == *((unsigned long long *)(&y));
	if (x != x && y != y) return true;
	return x == y;
}

void check(double x, double y)
{
	volatile double r1, r2;

	r1 = hwround::add_up(x, y);
	r2 = nohwround::add_up(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "add_up error\n"; std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::add_down(x, y);
	r2 = nohwround::add_down(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "add_down error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::sub_up(x, y);
	r2 = nohwround::sub_up(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "sub_up error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::sub_down(x, y);
	r2 = nohwround::sub_down(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "sub_down error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::mul_up(x, y);
	r2 = nohwround::mul_up(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "mul_up error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::mul_down(x, y);
	r2 = nohwround::mul_down(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "mul_down error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::div_up(x, y);
	r2 = nohwround::div_up(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "div_up error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::div_down(x, y);
	r2 = nohwround::div_down(x, y);
	if (!samedouble(r1, r2)) {
		std::cout << "div_down error\n";
		std::cout << x << "\n";
		std::cout << y << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::sqrt_up(x);
	r2 = nohwround::sqrt_up(x);
	if (!samedouble(r1, r2)) {
		std::cout << "sqrt_up error\n";
		std::cout << x << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}

	r1 = hwround::sqrt_down(x);
	r2 = nohwround::sqrt_down(x);
	if (!samedouble(r1, r2)) {
		std::cout << "sqrt_down error\n";
		std::cout << x << "\n";
		std::cout << r1 << "\n";
		std::cout << r2 << "\n";
	}
}

int main() {
	double x, y;
	int i, j;
	unsigned long long t;

	double specials[11] = {
		0., 
		-0.,
		std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity(),
		(std::numeric_limits<double>::max)(),
		-(std::numeric_limits<double>::max)(),
		(std::numeric_limits<double>::min)(),
		-(std::numeric_limits<double>::min)(),
		std::numeric_limits<double>::denorm_min(),
		-std::numeric_limits<double>::denorm_min()
	};
	specials[10] = specials[2] + specials[3]; // making NaN

	boost::variate_generator<boost::mt19937, boost::uniform_int<unsigned long long> > rand(boost::mt19937(time(0)), boost::uniform_int<unsigned long long>(0, -1));

	std::cout.precision(17);

	// cause overflow of intermediate variable in twoproduct
	// x = 6.929001713869936e+236;
	// y = 2.5944475251952003e+71;
	// check(x, y);

	// cause overflow of intermediate variable in twosum
	// x = 3.5630624444874539e+307;
	// y = -1.7976931348623157e+308;
	// check(x, y);

	// check general-general case

	for (i=0; i<NT; i++) {
		t = rand();
		x = *((double*)(&t));
		t = rand();
		y = *((double*)(&t));
		check(x, y);
	}

	// check general-special case

	for (i=0; i<NT; i++) {
		t = rand();
		x = *((double*)(&t));
		for (j=0; j<11; j++) {
			y = specials[j];
			check(x, y);
			check(y, x);
		}
	}

	// check special-special case

	for (i=0; i<11; i++) {
		x = specials[i];
		for (j=0; j<11; j++) {
			y = specials[j];
			check(x, y);
		}
	}
}
