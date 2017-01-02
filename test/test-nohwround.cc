/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

/*
 * test program for rounded math emulation
 *  compare add, sub, mul, div and sqrt under
 *  - changing rounding mode
 *  - emulating rounded math using twosum and twoproduct
 */

#include <fenv.h>
#include <boost/random.hpp>

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

	static void twosum(const double& a, const double& b, double& x, double& y) {
		double tmp;
		x = a + b;
		tmp = x - a;
		y = (a - (x - tmp)) + (b - tmp);
	}

	static void split(const double& a, double& x, double& y) {
		static const double sigma = (double)((1L << 27) + 1);
		double tmp;
		tmp = a * sigma;
		x = tmp - (tmp - a);
		y = a - x;
	}

	static void twoproduct(const double& a, const double& b, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);

		double na, nb, a1, a2, b1, b2;

		x = a * b;
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
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
		y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
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

		a = fabs(x);
		if (a >= th1) return x + a * c1;
		if (a < th2) return x + c2;
		c = c3 * x;
		e = c1 * fabs(c);
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

		a = fabs(x);
		if (a >= th1) return x - a * c1;
		if (a < th2) return x - c2;
		c = c3 * x;
		e = c1 * fabs(c);
		return (c - e) * c4;
	}


	static double add_up(const double& x, const double& y) {
		double r, r2;

		twosum(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return r;
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<double>::max)();
		}
		#ifdef NANCHECK
		if (r != r) return std::numeric_limits<double>::infinity();
		#endif

		if (r2 > 0.) {
			return succ(r);
		}

		return r;
	}

	static double add_down(const double& x, const double& y) {
		double r, r2;

		twosum(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<double>::max)();
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}
		#ifdef NANCHECK
		if (r != r) return -std::numeric_limits<double>::infinity();
		#endif

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
			return -(std::numeric_limits<double>::max)();
		}
		#ifdef NANCHECK
		if (r != r) return std::numeric_limits<double>::infinity();
		#endif

		if (r2 > 0.) {
			return succ(r);
		}

		return r;
	}

	static double sub_down(const double& x, const double& y) {
		double r, r2;

		twosum(x, -y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<double>::max)();
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}
		#ifdef NANCHECK
		if (r != r) return -std::numeric_limits<double>::infinity();
		#endif

		if (r2 < 0.) {
			return pred(r);
		}

		return r;
	}

	static double mul_up(const double& x, const double& y) {
		double r, r2;
		double x1, y1;
		int x2, y2;
		double s, s2, t;
		static const double th = ldexp(1., -970);

		if (x == 0. || y == 0.) return 0.;

		twoproduct(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return r;
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<double>::max)();
		}
		#ifdef NANCHECK
		if (r != r) return std::numeric_limits<double>::infinity();
		#endif

		if (fabs(r) >= th) {
			if (r2 > 0.) return succ(r);
			return r;
		} else {
			x1 = frexp(x, &x2);
			y1 = frexp(y, &y2);
			twoproduct(x1, y1, s, s2);
			t = ldexp(r, - x2 - y2);
			if ( t < s || (t == s && s2 > 0.)) {
				return succ(r);
			}
			return r;
		}
	}

	static double mul_down(const double& x, const double& y) {
		double r, r2;
		double x1, y1;
		int x2, y2;
		double s, s2, t;
		static const double th = ldexp(1., -970);

		if (x == 0. || y == 0.) return 0.;

		twoproduct(x, y, r, r2);
		if (r == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<double>::max)();
		} else if (r == -std::numeric_limits<double>::infinity()) {
			return r;
		}
		#ifdef NANCHECK
		if (r != r) return -std::numeric_limits<double>::infinity();
		#endif

		if (fabs(r) >= th) {
			if (r2 < 0.) return pred(r);
			return r;
		} else {
			x1 = frexp(x, &x2);
			y1 = frexp(y, &y2);
			twoproduct(x1, y1, s, s2);
			t = ldexp(r, - x2 - y2);
			if ( t > s || (t == s && s2 < 0.)) {
				return pred(r);
			}
			return r;
		}
	}

	static double div_up(const double& x, const double& y) {
		double r, r2;
		double xn, yn, d;
		static const double th1 = ldexp(1., -970);
		static const double th2 = ldexp(1., 919);
		static const double c1 = ldexp(1., 104);
		static const double c2 = ldexp(1., -1074);
		bool flag = false;

		if (x == 0. ) return x / y;
		if (x != x || y != y) return x / y;

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
				flag = true;
			}
		}

		d = xn / yn;

		if (d == std::numeric_limits<double>::infinity()) {
			return d;
		} else if (d == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<double>::max)();
		}
		#ifdef NANCHECK
		if (d != d) return std::numeric_limits<double>::infinity();
		#endif

		if (flag) {
			if (xn < 0.) return 0.;
			else return c2;
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
		static const double th1 = ldexp(1., -970);
		static const double th2 = ldexp(1., 919);
		static const double c1 = ldexp(1., 104);
		static const double c2 = ldexp(1., -1074);
		bool flag = false;

		if (x == 0. ) return x / y;
		if (x != x || y != y) return x / y;

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
				flag = true;
			}
		}

		d = xn / yn;

		if (d == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<double>::max)();
		} else if (d == -std::numeric_limits<double>::infinity()) {
			return d;
		}
		#ifdef NANCHECK
		if (d != d) return -std::numeric_limits<double>::infinity();
		#endif

		if (flag) {
			if (xn < 0.) return -c2;
			else return 0.;
		}

		twoproduct(d, yn, r, r2);
		if ( r > xn || ((r == xn) && r2 > 0.)) {
			return pred(d);
		}
		return d;
	}

	static double sqrt_up(const double& x) {
		double r, r2, d;
		static const double th1 = ldexp(1., -970);
		static const double c1 = ldexp(1., 104);
		static const double c2 = ldexp(1., 52);

		d = sqrt(x);

		if (x < th1) {
			double d2, x2;
			x2 = x * c1;
			d2 = d * c2;
			twoproduct(d2, d2, r, r2);
			if ( r < x2 || (r == x2) && r2 < 0.) {
				return succ(d);
			}
			return d;
		}

		twoproduct(d, d, r, r2);
		if ( r < x || (r == x) && r2 < 0.) {
			return succ(d);
		}
		return d;
	}

	static double sqrt_down(const double& x) {
		double r, r2, d;
		static const double th1 = ldexp(1., -970);
		static const double c1 = ldexp(1., 104);
		static const double c2 = ldexp(1., 52);

		d = sqrt(x);

		if (x < th1) {
			double d2, x2;
			x2 = x * c1;
			d2 = d * c2;
			twoproduct(d2, d2, r, r2);
			if ( r > x2 || (r == x2) && r2 > 0.) {
				return pred(d);
			}
			return d;
		}

		twoproduct(d, d, r, r2);
		if ( r > x || (r == x) && r2 > 0.) {
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

int main() {
	double x, y;
	volatile double r1, r2;
	int i;
	unsigned long long t;

	boost::variate_generator<boost::mt19937, boost::uniform_int<unsigned long long> > rand(boost::mt19937(time(0)), boost::uniform_int<unsigned long long>(0, -1));

	std::cout.precision(17);

	for (i=0; i<100000000; i++) {
		t = rand();
		x = *((double*)(&t));
		t = rand();
		y = *((double*)(&t));

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
}
