#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <fenv.h>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace kv {

template <class T> class rop {
	public:

	static T add_up(const T& x, const T& y) {
		return x + y;
	}

	static T add_down(const T& x, const T& y) {
		return x + y;
	}

	static T sub_up(const T& x, const T& y) {
		return x - y;
	}

	static T sub_down(const T& x, const T& y) {
		return x - y;
	}

	static T mul_up(const T& x, const T& y) {
		return x * y;
	}

	static T mul_down(const T& x, const T& y) {
		return x * y;
	}

	static T div_up(const T& x, const T& y) {
		return x / y;
	}

	static T div_down(const T& x, const T& y) {
		return x / y;
	}

	static T sqrt_up(const T& x) {
		return sqrt(x);
	}

	static T sqrt_down(const T& x) {
		return sqrt(x);
	}

	static void begin() {
	}

	static void finish() {
	}
};

template <class T> class interval {
	T inf;
	T sup;

	public:

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

		rop<T>::begin();
		r.inf = rop<T>::add_down(x.inf, y.inf);
		r.sup = rop<T>::add_up(x.sup, y.sup);
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator+(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::add_down(x.inf, y);
		r.sup = rop<T>::add_up(x.sup, y);
		rop<T>::finish();
		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator+(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::add_down(x, y.inf);
		r.sup = rop<T>::add_up(x, y.sup);
		rop<T>::finish();

		return r;
	}

	friend interval& operator+=(interval& x, const interval& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator+=(interval& x, const C& y) {

		rop<T>::begin();
		x.inf = rop<T>::add_down(x.inf + y);
		x.sup = rop<T>::add_up(x.sup + y);
		rop<T>::finish();

		return x;
	}

	friend interval operator-(const interval& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(x.inf, y.sup);
		r.sup = rop<T>::sub_up(x.sup, y.inf);
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator-(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(x.inf, y);
		r.sup = rop<T>::sub_up(x.sup, y);
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator-(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(x, y.sup);
		r.sup = rop<T>::sub_up(x, y.inf);
		rop<T>::finish();

		return r;
	}

	friend interval& operator-=(interval& x, const interval& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval& >::type operator-=(interval& x, const C& y) {

		rop<T>::begin();
		x.inf = rop<T>::sub_down(x.inf, y);
		x.sup = rop<T>::sub_up(x.sup, y);
		rop<T>::finish();

		return x;
	}

	friend interval operator-(const interval& x) {
		interval r;

		rop<T>::begin();
		r.sup = - x.inf;
		r.inf = - x.sup;
		rop<T>::finish();

		return r;
	}

	friend interval operator*(const interval& x, const interval& y) {
		interval r;
		T tmp;

		rop<T>::begin();
		if (x.inf >= 0.) {
			if (y.inf >= 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.inf);
				r.sup = rop<T>::mul_up(x.sup, y.sup);
			} else if (y.sup <= 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				r.sup = rop<T>::mul_up(x.inf, y.sup);
			} else {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				r.sup = rop<T>::mul_up(x.sup, y.sup);
			}
		} else if (x.sup <= 0.) {
			if (y.inf >= 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				r.sup = rop<T>::mul_up(x.sup, y.inf);
			} else if (y.sup <= 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.sup);
				r.sup = rop<T>::mul_up(x.inf, y.inf);
			} else {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				r.sup = rop<T>::mul_up(x.inf, y.inf);
			}
		} else {
			if (y.inf >= 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				r.sup = rop<T>::mul_up(x.sup, y.sup);
			} else if (y.sup <= 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				r.sup = rop<T>::mul_up(x.inf, y.inf);
			} else {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				tmp = rop<T>::mul_down(x.sup, y.inf);
				if (tmp < r.inf) r.inf = tmp;
				r.sup = rop<T>::mul_up(x.inf, y.inf);
				tmp = rop<T>::mul_up(x.sup, y.sup);
				if (tmp > r.sup) r.sup = tmp;
			}
		}
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator*(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		if (y >= 0.) {
			r.inf = rop<T>::mul_down(x.inf, y);
			r.sup = rop<T>::mul_up(x.sup, y);
		} else {
			r.inf = rop<T>::mul_down(x.sup, y);
			r.sup = rop<T>::mul_up(x.inf, y);
		}
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator*(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		if (x >= 0.) {
			r.inf = rop<T>::mul_down(x, y.inf);
			r.sup = rop<T>::mul_up(x, y.sup);
		} else {
			r.inf = rop<T>::mul_down(x, y.sup);
			r.sup = rop<T>::mul_up(x, y.inf);
		}
		rop<T>::finish();

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

		rop<T>::begin();
		if (y.inf > 0.) {
			if (x.inf >= 0.) {
				r.inf = rop<T>::div_down(x.inf, y.sup);
				r.sup = rop<T>::div_up(x.sup, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = rop<T>::div_down(x.inf, y.inf);
				r.sup = rop<T>::div_up(x.sup, y.sup);
			} else {
				r.inf = rop<T>::div_down(x.inf, y.inf);
				r.sup = rop<T>::div_up(x.sup, y.inf);
			}
		} else if (y.sup < 0.) {
			if (x.inf >= 0.) {
				r.inf = rop<T>::div_down(x.sup, y.sup);
				r.sup = rop<T>::div_up(x.inf, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = rop<T>::div_down(x.sup, y.inf);
				r.sup = rop<T>::div_up(x.inf, y.sup);
			} else {
				r.inf = rop<T>::div_down(x.sup, y.sup);
				r.sup = rop<T>::div_up(x.inf, y.sup);
			}
		} else {
			rop<T>::finish();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator/(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		if (y > 0.) {
			r.inf = rop<T>::div_down(x.inf, y);
			r.sup = rop<T>::div_up(x.sup, y);
		} else if (y < 0.) {
			r.inf = rop<T>::div_down(x.sup, y);
			r.sup = rop<T>::div_up(x.inf, y);
		} else {
			rop<T>::finish();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::finish();

		return r;
	}

	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, interval >::type operator/(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		if (y.inf > 0. || y.sup < 0.) {
			if (x >= 0.) {
				r.inf = rop<T>::div_down(x, y.sup);
				r.sup = rop<T>::div_up(x, y.inf);
			} else {
				r.inf = rop<T>::div_down(x, y.inf);
				r.sup = rop<T>::div_up(x, y.sup);
			}
		} else {
			rop<T>::finish();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::finish();

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
		T tmp1, tmp2;

		rop<T>::begin();
		tmp1 = rop<T>::sqrt_down(x.inf);
		tmp2 = rop<T>::sqrt_up(x.sup);
		rop<T>::finish();

		return interval(tmp1, tmp2);
	}

	const T& lower() const {
		return inf;
	}
	const T& upper() const {
		return sup;
	}
	T& lower() {
		return inf;
	}
	T& upper() {
		return sup;
	}

	static interval whole() {
		return interval(-std::numeric_limits<T>::infinity(),
		                 std::numeric_limits<T>::infinity() );
	}

	static interval hull(const T& x, const T& y) {
		if (x < y) return interval(x, y);
		else return interval(y, x);
	}

	friend T width(const interval& x) {
		T tmp;

		rop<T>::begin();
		tmp = rop<T>::sub_up(x.sup, x.inf);
		rop<T>::finish();

		return tmp;
	}

	friend T rad(const interval& x) {
		T tmp;

		rop<T>::begin();
		tmp = rop<T>::mul_up(rop<T>::sub_up(x.sup, x.inf), 0.5);
		rop<T>::finish();

		return tmp;
	}

	friend T median(const interval& x) {
		return (x.inf + x.sup) / 2.;
	}

	friend T mid(const interval& x) {
		return median(x);
	}

	friend T norm(const interval& x) {
		if (x.inf >= 0.) return x.sup;
		if (x.sup <= 0.) return -x.inf;
		T tmp = -x.inf;
		if (x.sup > tmp) tmp = x.sup;
		return tmp;
	}

	friend T mag(const interval& x) {
		return norm(x);
	}

	friend T mig(const interval& x) {
		if (zero_in(x)) return 0.;
		if (x.inf > 0.) return x.inf;
		else return -x.sup;
	}

	friend interval abs(const interval& x) {
		if (x.inf >= 0.) return x;
		if (x.sup <= 0.) return -x;
		T tmp = -x.inf;
		if (x.sup > tmp) tmp = x.sup;
		return interval(0., tmp);
	}

	friend bool in(const T& a, const interval& x) {
		return x.inf <= a && a <= x.sup;
	}

	friend bool zero_in(const interval& x) {
		return x.inf <= 0. && 0. <= x.sup;
	}

	friend bool subset(const interval& x, const interval& y) {
		return y.inf <= x.inf && x.sup <= y.sup;
	}

	friend bool proper_subset(const interval& x, const interval& y) {
		return y.inf < x.inf && x.sup < y.sup;
	}

	friend bool overlap(const interval& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y.inf > tmp1) tmp1 = y.inf;
		tmp2 = x.sup;
		if (y.sup < tmp2) tmp2 = y.sup;
		return tmp1 <= tmp2;
	}

	friend interval intersect(const interval& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y.inf > tmp1) tmp1 = y.inf;
		tmp2 = x.sup;
		if (y.sup < tmp2) tmp2 = y.sup;
		return interval(tmp1, tmp2);
	}

	friend interval hull(const interval& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y.inf < tmp1) tmp1 = y.inf;
		tmp2 = x.sup;
		if (y.sup > tmp2) tmp2 = y.sup;
		return interval(tmp1, tmp2);
	}

	friend interval hull(const interval& x, const T& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y < tmp1) tmp1 = y;
		tmp2 = x.sup;
		if (y > tmp2) tmp2 = y;
		return interval(tmp1, tmp2);
	}

	friend interval hull(const T& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x;
		if (y.inf < tmp1) tmp1 = y.inf;
		tmp2 = x;
		if (y.sup > tmp2) tmp2 = y.sup;
		return interval(tmp1, tmp2);
	}

	friend bool operator<(const interval& x, const interval& y) {
		return x.sup < y.inf;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator<(const interval& x, const C& y) {
		return x.sup < (T)y;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator<(const C& x, const interval& y) {
		return (T)x < y.inf;
	}

	friend bool operator<=(const interval& x, const interval& y) {
		return x.sup <= y.inf;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator<=(const interval& x, const C& y) {
		return x.sup <= (T)y;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator<=(const C& x, const interval& y) {
		return (T)x <= y.inf;
	}

	friend bool operator>(const interval& x, const interval& y) {
		return x.sup > y.inf;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator>(const interval& x, const C& y) {
		return x.sup > (T)y;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator>(const C& x, const interval& y) {
		return (T)x > y.inf;
	}

	friend bool operator>=(const interval& x, const interval& y) {
		return x.sup >= y.inf;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator>=(const interval& x, const C& y) {
		return x.sup >= (T)y;
	}
	template <class C> friend typename boost::enable_if< boost::is_convertible<C, T>, bool>::type operator>=(const C& x, const interval& y) {
		return (T)x >= y.inf;
	}

	friend interval division_part1(const interval& x, const interval& y, bool& parted) {
		interval r;

		parted = false;

		if (y.inf > 0. || y.sup < 0.) {
			return x / y;
		}

		if (y.inf == 0. || y.sup == 0.) {
			throw std::range_error("interval: division by 0");
		}

		rop<T>::begin();
		if (y.inf < 0. && y.sup == 0.) {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf >= 0.) {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.inf, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = rop<T>::div_down(x.sup, y.inf);
				r.sup = std::numeric_limits<T>::infinity();
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		} else if (y.inf == 0. && y.sup > 0.) {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf >= 0.) {
				r.inf = rop<T>::div_down(x.inf, y.sup);
				r.sup = std::numeric_limits<T>::infinity();
			} else if (x.sup <= 0.) {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.sup, y.sup);
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		} else {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf > 0.) {
				parted = true;
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.inf, y.inf);
			} else if (x.sup < 0.) {
				parted = true;
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.sup, y.sup);
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		}
		rop<T>::finish();

		return r;
	}

	friend interval division_part2(const interval& x, const interval& y, bool parted = true) {
		interval r;

		if (y.inf > 0. || y.sup < 0.) {
			return x / y;
		}

		if (y.inf == 0. || y.sup == 0.) {
			throw std::range_error("interval: division by 0");
		}

		rop<T>::begin();
		if (y.inf < 0. && y.sup == 0.) {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf >= 0.) {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.inf, y.inf);
			} else if (x.sup <= 0.) {
				r.inf = rop<T>::div_down(x.sup, y.inf);
				r.sup = std::numeric_limits<T>::infinity();
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		} else if (y.inf == 0. && y.sup > 0.) {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf >= 0.) {
				r.inf = rop<T>::div_down(x.inf, y.sup);
				r.sup = std::numeric_limits<T>::infinity();
			} else if (x.sup <= 0.) {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = rop<T>::div_up(x.sup, y.sup);
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		} else {
			if (x.inf == 0. && x.sup == 0) {
				r.inf = 0.;
				r.sup = 0.;
			} else if (x.inf > 0.) {
				r.inf = rop<T>::div_down(x.inf, y.sup);
				r.sup = std::numeric_limits<T>::infinity();
			} else if (x.sup < 0.) {
				r.inf = rop<T>::div_down(x.sup, y.inf);
				r.sup = std::numeric_limits<T>::infinity();
			} else {
				r.inf = -std::numeric_limits<T>::infinity();
				r.sup = std::numeric_limits<T>::infinity();
			}
		}
		rop<T>::finish();

		return r;
	}

};

} // namespace kv

#endif // INTERVAL_HPP
