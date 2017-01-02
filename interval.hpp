#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fenv.h>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace kv {

	template <class T> inline T add_up(const T& x, const T& y) {
		return x + y;
	}

	template <class T> inline T add_down(const T& x, const T& y) {
		return x + y;
	}

	template <class T> inline T sub_up(const T& x, const T& y) {
		return x - y;
	}

	template <class T> inline T sub_down(const T& x, const T& y) {
		return x - y;
	}

	template <class T> inline T mul_up(const T& x, const T& y) {
		return x * y;
	}

	template <class T> inline T mul_down(const T& x, const T& y) {
		return x * y;
	}

	template <class T> inline T div_up(const T& x, const T& y) {
		return x / y;
	}

	template <class T> inline T div_down(const T& x, const T& y) {
		return x / y;
	}

	template <class T> inline T sqrt_up(const T& x) {
		return sqrt(x);
	}

	template <class T> inline T sqrt_down(const T& x) {
		return sqrt(x);
	}

	inline void roundnear() {
		fesetround(FE_TONEAREST);
	}

	inline void rounddown() {
		fesetround(FE_DOWNWARD);
	}

	inline void roundup() {
		fesetround(FE_UPWARD);
	}

	inline void roundchop() {
		fesetround(FE_TOWARDZERO);
	}

	template <> inline double add_up<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 + y1;
		roundnear();
		return r;
	}

	template <> inline double add_down<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 + y1;
		roundnear();
		return r;
	}

	template <> inline double sub_up<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 - y1;
		roundnear();
		return r;
	}

	template <> inline double sub_down<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 - y1;
		roundnear();
		return r;
	}

	template <> inline double mul_up<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 * y1;
		roundnear();
		return r;
	}

	template <> inline double mul_down<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 * y1;
		roundnear();
		return r;
	}

	template <> inline double div_up<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		roundup();
		r = x1 / y1;
		roundnear();
		return r;
	}

	template <> inline double div_down<double>(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		rounddown();
		r = x1 / y1;
		roundnear();
		return r;
	}

	template <> inline double sqrt_up<double>(const double& x) {
		volatile double r, x1 = x;
		roundup();
		r = sqrt(x1);
		roundnear();
		return r;
	}

	template <> inline double sqrt_down<double>(const double& x) {
		volatile double r, x1 = x;
		rounddown();
		r = sqrt(x1);
		roundnear();
		return r;
	}

template <class T> class interval {
	public:
	T inf;
	T sup;

	interval() {
	}

	template <class C> interval(const C& x, typename boost::enable_if< boost::is_convertible<C, T> >::type* =0) {
		inf = x;
		sup = x;
	}

	template <class C1, class C2> interval(const C1& x, const C2& y, typename boost::enable_if_c< boost::is_convertible<C1, T>::value && boost::is_convertible<C2, T>::value >::type* =0) {
		inf = x;
		sup = y;
	}

	friend interval operator+(const interval& x, const interval& y) {
		interval r;

		r.inf = add_down(x.inf, y.inf);
		r.sup = add_up(x.sup, y.sup);

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator+(const interval& x, const C& y) {
		interval r;

		r.inf = add_down(x.inf, y);
		r.sup = add_up(x.sup, y);

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator+(const C& x, const interval& y) {
		interval r;

		r.inf = add_down(x, y.inf);
		r.sup = add_up(x, y.sup);

		return r;
	}

	friend interval& operator+=(interval& x, const interval& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator+=(interval& x, const C& y) {
		x.inf = add_down(x.inf + y);
		x.sup = add_up(x.sup + y);
		return x;
	}

	friend interval operator-(const interval& x, const interval& y) {
		interval r;

		r.inf = sub_down(x.inf, y.sup);
		r.sup = sub_up(x.sup, y.inf);

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator-(const interval& x, const C& y) {
		interval r;

		r.inf = sub_down(x.inf, y);
		r.sup = sub_up(x.sup, y);

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator-(const C& x, const interval& y) {
		interval r;

		r.inf = sub_down(x, y.sup);
		r.sup = sub_up(x, y.inf);

		return r;
	}

	friend interval& operator-=(interval& x, const interval& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator-=(interval& x, const C& y) {
		x.inf = sub_down(x.inf, y);
		x.sup = sub_up(x.sup, y);
		return x;
	}

	friend interval operator-(const interval& x) {
		interval r;

		r.sup = - x.inf;
		r.inf = - x.sup;

		return r;
	}

	friend interval operator*(const interval& x, const interval& y) {
		interval r;
		T tmp;

		if (x.inf >= 0.) {
			if (y.inf >= 0.) {
				r.inf = mul_down(x.inf, y.inf);
				r.sup = mul_up(x.sup, y.sup);
			} else if (y.sup <= 0.) {
				r.inf = mul_down(x.sup, y.inf);
				r.sup = mul_up(x.inf, y.sup);
			} else {
				r.inf = mul_down(x.sup, y.inf);
				r.sup = mul_up(x.sup, y.sup);
			}
		} else if (x.sup <= 0.) {
			if (y.inf >= 0.) {
				r.inf = mul_down(x.inf, y.sup);
				r.sup = mul_up(x.sup, y.inf);
			} else if (y.sup <= 0.) {
				r.inf = mul_down(x.sup, y.sup);
				r.sup = mul_up(x.inf, y.inf);
			} else {
				r.inf = mul_down(x.inf, y.sup);
				r.sup = mul_up(x.inf, y.inf);
			}
		} else {
			if (y.inf >= 0.) {
				r.inf = mul_down(x.inf, y.sup);
				r.sup = mul_up(x.sup, y.sup);
			} else if (y.sup <= 0.) {
				r.inf = mul_down(x.sup, y.inf);
				r.sup = mul_up(x.inf, y.inf);
			} else {
				r.inf = mul_down(x.inf, y.sup);
				tmp = mul_down(x.sup, y.inf);
				if (tmp < r.inf) r.inf = tmp;
				r.sup = mul_up(x.inf, y.inf);
				tmp = mul_up(x.sup, y.sup);
				if (tmp > r.inf) r.inf = tmp;
			}
		}

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator*(const interval& x, const C& y) {
		interval r;

		if (y >= 0.) {
			r.inf = mul_down(x.inf, y);
			r.sup = mul_up(x.sup, y);
		} else {
			r.inf = mul_down(x.sup, y);
			r.sup = mul_up(x.inf, y);
		}

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator*(const C& x, const interval& y) {
		interval r;

		if (x >= 0.) {
			r.inf = mul_down(x, y.inf);
			r.sup = mul_up(x, y.sup);
		} else {
			r.inf = mul_down(x, y.sup);
			r.sup = mul_up(x, y.inf);
		}

		return r;
	}

	friend interval& operator*=(interval& x, const interval& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator*=(interval& x, const C& y) {
		x = x * y;
		return x;
	}

	friend interval operator/(const interval& x, const interval& y) {
		interval r;

		if (y.inf > 0.) {
			if (x.inf >= 0.) {
				r.inf = div_down(x.inf, y.sup);
				r.sup = div_up(x.sup, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = div_down(x.inf, y.inf);
				r.sup = div_up(x.sup, y.sup);
			} else {
				r.inf = div_down(x.inf, y.inf);
				r.sup = div_up(x.sup, y.inf);
			}
		} else if (y.sup < 0.) {
			if (x.inf >= 0.) {
				r.inf = div_down(x.sup, y.sup);
				r.sup = div_up(x.inf, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = div_down(x.sup, y.inf);
				r.sup = div_up(x.inf, y.sup);
			} else {
				r.inf = div_down(x.sup, y.sup);
				r.sup = div_up(x.inf, y.sup);
			}
		} else {
			throw std::range_error("interval: division by 0");
		}

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator/(const interval& x, const C& y) {
		interval r;

		if (y > 0.) {
			r.inf = div_down(x.inf, y);
			r.sup = div_up(x.sup, y);
		} else if (y < 0.) {
			r.inf = div_down(x.sup, y);
			r.sup = div_up(x.inf, y);
		} else {
			throw std::range_error("interval: division by 0");
		}

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator/(const C& x, const interval& y) {
		interval r;

		if (y.inf > 0. || y.sup < 0.) {
			if (x >= 0.) {
				r.inf = div_down(x, y.sup);
				r.sup = div_up(x, y.inf);
			} else {
				r.inf = div_down(x, y.inf);
				r.sup = div_up(x, y.sup);
			}
		} else {
			throw std::range_error("interval: division by 0");
		}

		return r;
	}

	friend interval& operator/=(interval& x, const interval& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator/=(interval& x, const C& y) {
		x = x / y;
		return x;
	}

	friend std::ostream& operator<<(std::ostream& s, const interval& x) {
		s << '[' << x.inf << ',' << x.sup << "]";
		return s;
	}

	friend interval sqrt(const interval& x) {
		interval r(sqrt_down(x.inf), sqrt_up(x.sup));
		return r;
	}
};

} // namespace kv

#endif // INTERVAL_HPP
