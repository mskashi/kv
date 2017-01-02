#ifndef RDOUBLE_HPP
#define RDOUBLE_HPP

namespace kv {

template <> class rop <double> {
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

	static void begin() {
	}

	static void finish() {
	}
};

} // namespace kv

#endif // RDOUBLE_HPP
