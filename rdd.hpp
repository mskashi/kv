#ifndef RDD_HPP
#define RDD_HPP

namespace kv {

template <> class rop <dd> {
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

	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		double a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;
		x = a * b;
		dd::split(a, a1, a2);
		dd::split(b, b1, b2);
		roundup();
		v1 = x; v3 = a1; v4 = a2; v5 =b1; v6 = b2; 
		v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) - v4 * v6;
		y = v2;
		roundnear();
	}

	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		double a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;
		x = a * b;
		dd::split(a, a1, a2);
		dd::split(b, b1, b2);
		rounddown();
		v1 = x; v3 = a1; v4 = a2; v5 =b1; v6 = b2; 
		v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) - v4 * v6;
		y = v2;
		roundnear();
	}

	static dd add_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, y.a1, z1, z2);
		roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, y.a1, z1, z2);
		rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);
		roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);
		rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_up(x.a1, y.a1, z1, z2);
		roundup();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_down(x.a1, y.a1, z1, z2);
		rounddown();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

		if (y >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			rounddown();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			roundup();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		} else {
			twoproduct_down(-z1, y.a1, z3, z4);
			roundup();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			rounddown();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		}
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd div_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

		if (y >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			roundup();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			rounddown();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		} else {
			twoproduct_up(-z1, y.a1, z3, z4);
			rounddown();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			roundup();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		}
		roundnear();
		dd::twosum(z1, z2, z3, z4);
		z1 = z3; z2 = z4;

		return dd(z1, z2);
	}

	static dd sqrt_up(const dd& x) {
		dd r, r2;

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_up(x, r);

		if (r > r2) return r;
		else return r2;
	}

	static dd sqrt_down(const dd& x) {
		dd r, r2;

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_down(x, r);

		if (r > r2) return r2;
		else return r;
	}

	static void begin() {
	}

	static void finish() {
	}

	static void print_up(const dd& x, std::ostream& s) {
		char format;
		if (s.flags() & s.scientific) {
			if (s.flags() & s.fixed) {
				format = 'g';
			} else {
				format = 'e';
			}
		} else {
			if (s.flags() & s.fixed) {
				format = 'f';
			} else {
				format = 'g';
			}
		}
		s << dd::ddtostring(x, s.precision(), format, 1);
	}

	static void print_down(const dd& x, std::ostream& s) {
		char format;
		if (s.flags() & s.scientific) {
			if (s.flags() & s.fixed) {
				format = 'g';
			} else {
				format = 'e';
			}
		} else {
			if (s.flags() & s.fixed) {
				format = 'f';
			} else {
				format = 'g';
			}
		}
		s << dd::ddtostring(x, s.precision(), format, -1);
	}

	static dd fromstring_up(const std::string& s) {
		return dd::stringtodd(s, 1);
	}

	static dd fromstring_down(const std::string& s) {
		return dd::stringtodd(s, -1);
	}
};

} // namespace kv

#endif // RDD_HPP
