/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef INTERVAL_HPP
#define INTERVAL_HPP

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <string>
#include <sstream>

#include <kv/convert.hpp>
#include <kv/constants.hpp>


namespace kv {

template <class T> struct rop {

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

	static void end() {
	}

	static void print_up(const T& x, std::ostream& s) {
		s << x;
	}

	static void print_down(const T& x, std::ostream& s) {
		s << x;
	}

	static T fromstring_up(const std::string& s) {
		std::istringstream is(s);
		T r;
		is >> r;
		return r;
	}

	static T fromstring_down(const std::string& s) {
		std::istringstream is(s);
		T r;
		is >> r;
		return r;
	}
};


template <class T> class interval;
template <class C, class T> struct convertible<C, interval<T> > {
	static const bool value = convertible<C, T>::value || boost::is_same<C, interval<T> >::value || boost::is_convertible<C, std::string>::value;
};
template <class C, class T> struct acceptable_n<C, interval<T> > {
	static const bool value = convertible<C, T>::value && (! boost::is_convertible<C, std::string>::value);
};
template <class C, class T> struct acceptable_s<C, interval<T> > {
	static const bool value = boost::is_convertible<C, std::string>::value;
};


template <class T> class interval {
	T inf;
	T sup;

	public:

	typedef T base_type;

	interval() {
		inf = 0.;
		sup = 0.;
	}

	template <class C> explicit interval(const C& x, typename boost::enable_if_c< acceptable_n<C, interval>::value >::type* =0) {
		inf = x;
		sup = x;
	}

	template <class C1, class C2> interval(const C1& x, const C2& y, typename boost::enable_if_c< acceptable_n<C1, interval>::value && acceptable_n<C2, interval>::value >::type* =0) {
		inf = x;
		sup = y;
	}

	template <class C> explicit interval(const C& x, typename boost::enable_if_c< acceptable_s<C, interval>::value >::type* =0) {
		inf = rop<T>::fromstring_down(x);
		sup = rop<T>::fromstring_up(x);
	}

	template <class C1, class C2> interval(const C1& x, const C2& y, typename boost::enable_if_c< acceptable_s<C1, interval>::value && acceptable_s<C2, interval>::value >::type* =0) {
		inf = rop<T>::fromstring_down(x);
		sup = rop<T>::fromstring_up(y);
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, interval>::value, interval& >::type operator=(const C& x) {
		inf = x;
		sup = x;
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_s<C, interval>::value, interval& >::type operator=(const C& x) {
		inf = rop<T>::fromstring_down(x);
		sup = rop<T>::fromstring_up(x);
		return *this;
	}

	friend interval operator+(const interval& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::add_down(x.inf, y.inf);
		r.sup = rop<T>::add_up(x.sup, y.sup);
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator+(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::add_down(x.inf, T(y));
		r.sup = rop<T>::add_up(x.sup, T(y));
		rop<T>::end();
		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval>::type operator+(const interval& x, const C& y) {
		return x + interval(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator+(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::add_down(T(x), y.inf);
		r.sup = rop<T>::add_up(T(x), y.sup);
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval>::type operator+(const C& x, const interval& y) {
		return interval(x) + y;
	}

	friend interval& operator+=(interval& x, const interval& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval& >::type operator+=(interval& x, const C& y) {

		rop<T>::begin();
		x.inf = rop<T>::add_down(x.inf, T(y));
		x.sup = rop<T>::add_up(x.sup, T(y));
		rop<T>::end();

		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval& >::type operator+=(interval& x, const C& y) {
		x = x + interval(y);
		return x;
	}

	friend interval operator-(const interval& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(x.inf, y.sup);
		r.sup = rop<T>::sub_up(x.sup, y.inf);
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator-(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(x.inf, T(y));
		r.sup = rop<T>::sub_up(x.sup, T(y));
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator-(const interval& x, const C& y) {
		return x - interval(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator-(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		r.inf = rop<T>::sub_down(T(x), y.sup);
		r.sup = rop<T>::sub_up(T(x), y.inf);
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator-(const C& x, const interval& y) {
		return interval(x) - y;
	}

	friend interval& operator-=(interval& x, const interval& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval& >::type operator-=(interval& x, const C& y) {

		rop<T>::begin();
		x.inf = rop<T>::sub_down(x.inf, T(y));
		x.sup = rop<T>::sub_up(x.sup, T(y));
		rop<T>::end();

		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval& >::type operator-=(interval& x, const C& y) {
		x = x - interval(y);
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

		rop<T>::begin();
		if (x.inf > 0.) {
			if (y.inf > 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.inf);
				r.sup = rop<T>::mul_up(x.sup, y.sup);
			} else if (y.sup < 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				r.sup = rop<T>::mul_up(x.inf, y.sup);
			} else {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				if (r.inf != r.inf) r = 0.;
				r.sup = rop<T>::mul_up(x.sup, y.sup);
				if (r.sup != r.sup) r = 0.;
			}
		} else if (x.sup < 0.) {
			if (y.inf > 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				r.sup = rop<T>::mul_up(x.sup, y.inf);
			} else if (y.sup < 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.sup);
				r.sup = rop<T>::mul_up(x.inf, y.inf);
			} else {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				if (r.inf != r.inf) r = 0.;
				r.sup = rop<T>::mul_up(x.inf, y.inf);
				if (r.sup != r.sup) r = 0.;
			}
		} else {
			if (y.inf > 0.) {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				if (r.inf != r.inf) r = 0.;
				r.sup = rop<T>::mul_up(x.sup, y.sup);
				if (r.sup != r.sup) r = 0.;
			} else if (y.sup < 0.) {
				r.inf = rop<T>::mul_down(x.sup, y.inf);
				if (r.inf != r.inf) r = 0.;
				r.sup = rop<T>::mul_up(x.inf, y.inf);
				if (r.sup != r.sup) r = 0.;
			} else {
				r.inf = rop<T>::mul_down(x.inf, y.sup);
				if (r.inf != r.inf) r = 0.;
				tmp = rop<T>::mul_down(x.sup, y.inf);
				if (tmp != tmp) r = 0.;
				if (tmp < r.inf) r.inf = tmp;
				r.sup = rop<T>::mul_up(x.inf, y.inf);
				if (r.sup != r.sup) r = 0.;
				tmp = rop<T>::mul_up(x.sup, y.sup);
				if (tmp != tmp) r = 0.;
				if (tmp > r.sup) r.sup = tmp;
			}
		}
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator*(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		if (y > 0.) {
			r.inf = rop<T>::mul_down(x.inf, T(y));
			r.sup = rop<T>::mul_up(x.sup, T(y));
		} else {
			r.inf = rop<T>::mul_down(x.sup, T(y));
			r.sup = rop<T>::mul_up(x.inf, T(y));
		}
		rop<T>::end();
		if (r.inf != r.inf) r = 0.;
		if (r.sup != r.sup) r = 0.;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator*(const interval& x, const C& y) {
		return x * interval(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator*(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		if (x >= 0.) {
			r.inf = rop<T>::mul_down(T(x), y.inf);
			r.sup = rop<T>::mul_up(T(x), y.sup);
		} else {
			r.inf = rop<T>::mul_down(T(x), y.sup);
			r.sup = rop<T>::mul_up(T(x), y.inf);
		}
		rop<T>::end();
		if (r.inf != r.inf) r = 0.;
		if (r.sup != r.sup) r = 0.;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator*(const C& x, const interval& y) {
		return interval(x) * y;
	}

	friend interval& operator*=(interval& x, const interval& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval& >::type operator*=(interval& x, const C& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval& >::type operator*=(interval& x, const C& y) {
		x = x * interval(y);
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
			rop<T>::end();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator/(const interval& x, const C& y) {
		interval r;

		rop<T>::begin();
		if (y > 0.) {
			r.inf = rop<T>::div_down(x.inf, T(y));
			r.sup = rop<T>::div_up(x.sup, T(y));
		} else if (y < 0.) {
			r.inf = rop<T>::div_down(x.sup, T(y));
			r.sup = rop<T>::div_up(x.inf, T(y));
		} else {
			rop<T>::end();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator/(const interval& x, const C& y) {
		return x / interval(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type operator/(const C& x, const interval& y) {
		interval r;

		rop<T>::begin();
		if (y.inf > 0. || y.sup < 0.) {
			if (x >= 0.) {
				r.inf = rop<T>::div_down(T(x), y.sup);
				r.sup = rop<T>::div_up(T(x), y.inf);
			} else {
				r.inf = rop<T>::div_down(T(x), y.inf);
				r.sup = rop<T>::div_up(T(x), y.sup);
			}
		} else {
			rop<T>::end();
			throw std::range_error("interval: division by 0");
		}
		rop<T>::end();

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval >::type operator/(const C& x, const interval& y) {
		return interval(x) / y;
	}

	friend interval& operator/=(interval& x, const interval& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, interval& >::type operator/=(interval& x, const C& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, interval& >::type operator/=(interval& x, const C& y) {
		x = x / interval(y);
		return x;
	}

	friend std::ostream& operator<<(std::ostream& s, const interval& x) {
		s << '[';
		rop<T>::print_down(x.inf, s);
		s << ',';
		rop<T>::print_up(x.sup, s);
		s << ']';
		return s;
	}

	friend interval sqrt(const interval& x) {
		T tmp1, tmp2;

		rop<T>::begin();
		tmp1 = rop<T>::sqrt_down(x.inf);
		tmp2 = rop<T>::sqrt_up(x.sup);
		rop<T>::end();

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

	template <class C1, class C2> typename boost::enable_if_c< acceptable_n<C1, interval>::value && acceptable_n<C2, interval>::value, void >::type assign(const C1& x, const C2& y) {
		inf = x;
		sup = y;
	}

	static interval whole() {
		return interval(-std::numeric_limits<T>::infinity(),
		                 std::numeric_limits<T>::infinity() );
	}

	template <class C1, class C2> static typename boost::enable_if_c< acceptable_n<C1, interval>::value && acceptable_n<C2, interval>::value, interval >::type hull(const C1& x, const C2& y) {
		if (x < y) return interval(x, y);
		else return interval(y, x);
	}

	static interval hull(const interval& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y.inf < tmp1) tmp1 = y.inf;
		tmp2 = x.sup;
		if (y.sup > tmp2) tmp2 = y.sup;
		return interval(tmp1, tmp2);
	}

	template <class C> static typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type hull(const interval& x, const C& y) {
		T tmp1, tmp2;
		tmp1 = x.inf;
		if (y < tmp1) tmp1 = y;
		tmp2 = x.sup;
		if (y > tmp2) tmp2 = y;
		return interval(tmp1, tmp2);
	}

	template <class C> static typename boost::enable_if_c< acceptable_n<C, interval>::value, interval >::type hull(const C& x, const interval& y) {
		T tmp1, tmp2;
		tmp1 = x;
		if (y.inf < tmp1) tmp1 = y.inf;
		tmp2 = x;
		if (y.sup > tmp2) tmp2 = y.sup;
		return interval(tmp1, tmp2);
	}


	friend T width(const interval& x) {
		T tmp;

		rop<T>::begin();
		tmp = rop<T>::sub_up(x.sup, x.inf);
		rop<T>::end();

		return tmp;
	}

	friend T rad(const interval& x) {
		T tmp;

		rop<T>::begin();
		tmp = rop<T>::mul_up(rop<T>::sub_up(x.sup, x.inf), T(0.5));
		rop<T>::end();

		return tmp;
	}

	friend T median(const interval& x) {
		return (x.inf + x.sup) * 0.5;
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
		if (zero_in(x)) return T(0.);
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

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool >::type in(const C& a, const interval& x) {
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


	friend bool operator<(const interval& x, const interval& y) {
		return x.sup < y.inf;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator<(const interval& x, const C& y) {
		return x.sup < T(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator<(const interval& x, const C& y) {
		return x.sup < rop<T>::fromstring_down(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator<(const C& x, const interval& y) {
		return T(x) < y.inf;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator<(const C& x, const interval& y) {
		return rop<T>::fromstring_up(y) < y.inf;
	}


	friend bool operator<=(const interval& x, const interval& y) {
		return x.sup <= y.inf;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator<=(const interval& x, const C& y) {
		return x.sup <= T(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator<=(const interval& x, const C& y) {
		return x.sup <= rop<T>::fromstring_down(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator<=(const C& x, const interval& y) {
		return T(x) <= y.inf;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator<=(const C& x, const interval& y) {
		return rop<T>::fromstring_up(x) <= y.inf;
	}


	friend bool operator>(const interval& x, const interval& y) {
		return x.inf > y.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator>(const interval& x, const C& y) {
		return x.inf > T(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator>(const interval& x, const C& y) {
		return x.inf > rop<T>::fromstring_up(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator>(const C& x, const interval& y) {
		return T(x) > y.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator>(const C& x, const interval& y) {
		return rop<T>::fromstring_down(x) > y.sup;
	}


	friend bool operator>=(const interval& x, const interval& y) {
		return x.inf >= y.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator>=(const interval& x, const C& y) {
		return x.inf >= T(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator>=(const interval& x, const C& y) {
		return x.inf >= rop<T>::fromstring_up(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator>=(const C& x, const interval& y) {
		return T(x) >= y.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator>=(const C& x, const interval& y) {
		return rop<T>::fromstring_down(x) >= y.sup;
	}


	friend bool operator==(const interval& x, const interval& y) {
		return x.inf == x.sup && x.sup == y.inf && y.inf == y.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator==(const interval& x, const C& y) {
		return x.inf == x.sup && x.sup == T(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator==(const interval& x, const C& y) {
		interval ytmp(y);
		return x.inf == x.sup && x.sup == ytmp.inf && ytmp.inf == ytmp.sup;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator==(const C& x, const interval& y) {
		return y.inf == y.sup && y.sup == T(x);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator==(const C& x, const interval& y) {
		interval xtmp(x);
		return xtmp.inf == xtmp.sup && xtmp.sup == y.inf && y.inf == y.sup;
	}


	friend bool operator!=(const interval& x, const interval& y) {
		return !overlap(x, y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator!=(const interval& x, const C& y) {
		return !overlap(x, interval(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator!=(const interval& x, const C& y) {
		return !overlap(x, interval(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, interval>::value, bool>::type operator!=(const C& x, const interval& y) {
		return !overlap(interval(x), y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, interval>::value, bool>::type operator!=(const C& x, const interval& y) {
		return !overlap(interval(x), y);
	}


	friend interval division_part1(const interval& x, const interval& y, bool& parted) {
		interval r;

		parted = false;

		if (y.inf > 0. || y.sup < 0.) {
			return x / y;
		}

		if (y.inf == 0. && y.sup == 0.) {
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
		rop<T>::end();

		return r;
	}

	friend interval division_part2(const interval& x, const interval& y, bool parted = true) {
		interval r;

		if (y.inf > 0. || y.sup < 0.) {
			return x / y;
		}

		if (y.inf == 0. && y.sup == 0.) {
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
		rop<T>::end();

		return r;
	}

	static interval pow_point(const T& x, int y) {
		interval r, xp;

		r = 1.;
		xp = (interval)x;

		while (y != 0) {
			if (y % 2 != 0) {
				r *= xp;
			}
			y /= 2;
			xp *= xp;
		}

		return r;
	}

	friend interval pow(const interval& x, int y) {
		interval r, xp;
		int a;

		if (y == 0) return interval(1.);

		a = (y >= 0) ? y : -y;

		if (a % 2 == 0 && in(0., x)) {
			r = hull(0., pow_point(mag(x), a));
		} else {
			r = hull(pow_point(x.lower(), a), pow_point(x.upper(), a));
		}

		if (y < 0) {
			r = 1. / r;
		}

		return r;
	}

	friend interval pow(const interval& x, const interval& y) {
		return exp(y * log(x));
	}

	// power by integer (i is passed by double)
	static interval ipower(const interval& x, double i) {
		double tmp;
		interval xp = x;
		interval r(1.);

		while (i != 0.) {
			i *= 0.5;
			tmp = floor(i);
			if (tmp != i) {
				i = tmp;
				r *= xp;
			}
			xp = xp * xp;
		}

		return r;
	}

	static interval exp_point(const T& x) {
		T x_i, x_f, tmp;
		interval r, y, remainder;
		int i;

		if (x == std::numeric_limits<T>::infinity()) {
			return interval((std::numeric_limits<T>::max()), std::numeric_limits<T>::infinity());
		}
		if (x == -std::numeric_limits<T>::infinity()) {
			return interval(0.);
		}

		remainder = interval((1./sqrt(constants<interval>::e())).lower(), sqrt(constants<interval>::e()).upper());

		using std::floor;
		if (x >= 0.) {
			x_i = floor(x);
			x_f = x - x_i;
			if (x_f >= 0.5) {
				x_f -= 1.;
				x_i += 1.;
			}
		} else {
			x_i = -floor(-x);
			x_f = x - x_i;
			if (x_f <= -0.5) {
				x_f += 1.;
				x_i -= 1.;
			}
		}

		r = 1.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= x_f;
			y /= (T)i;
			if (mag(y) * remainder.upper() < std::numeric_limits<T>::epsilon()) {
				r += y * remainder;
				break;
			} else {
				r += y;
			}
		}

		if (x_i >= 0.) {
			// r *= pow(constants<interval>::e(), (int)x_i);
			r *= ipower(constants<interval>::e(), (double)x_i);
		} else {
			// r /= pow(constants<interval>::e(), -(int)x_i);
			r /= ipower(constants<interval>::e(), -(double)x_i);
		}

		return r;
	}

	friend interval exp(const interval& I) {
		return interval(exp_point(I.lower()).lower(), exp_point(I.upper()).upper());
	}

	static interval expm1_origin(const T& x) {
		interval r, y, remainder;
		int i;

		remainder = interval((1./sqrt(constants<interval>::e())).lower(), sqrt(constants<interval>::e()).upper());

		r = 0.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= x;
			y /= (T)i;
			if (mag(y) * remainder.upper() < std::numeric_limits<T>::epsilon()) {
				r += y * remainder;
				break;
			} else {
				r += y;
			}
		}

		return r;
	}

	static interval expm1_point(const T& x) {
		if (x >= -0.5 && x <= 0.5) {
			return expm1_origin(x);
		} else {
			return exp_point(x) - 1.;
		}
	}

	friend interval expm1(const interval& I) {
		return interval(expm1_point(I.lower()).lower(), expm1_point(I.upper()).upper());
	}

	static T log_point(const T& x, int round) {
		interval tmp;
		T x2, x2m1;
		interval cinv;
		interval r;
		interval xn, xn2;
		interval sqrt2 = sqrt(interval((T)2.));
		int p_i;
		T p;
		int i;

		if (x == std::numeric_limits<T>::infinity()) {
			if (round == 1) {
				return std::numeric_limits<T>::infinity();
			} else {
				return (std::numeric_limits<T>::max)();
			}
		}
		if (x == 0.) {
			if (round == 1) {
				return -(std::numeric_limits<T>::max)();
			} else {
				return -std::numeric_limits<T>::infinity();
			}
		}

		using std::frexp;
		x2 = frexp(x, &p_i);
		p = (T)p_i;

		while (x2 > 4. * std::sqrt(2.) - 4.) {
			x2 *= 0.5;
			p += 1.;
		}
		while (x2 > 4. - 2. * std::sqrt(2.)) {
			tmp = x2 / sqrt2;
			if (round == -1) x2 = tmp.lower();
			else x2 = tmp.upper();
			p += 0.5;
		}
		while (x2 < 2. - std::sqrt(2.)) {
			x2 *= 2.;
			p -= 1.;
		}
		while (x2 < 2. * std::sqrt(2.) - 2.) {
			tmp = x2 * sqrt2;
			if (round == -1) x2 = tmp.lower();
			else x2 = tmp.upper();
			p -= 0.5;
		}

		x2m1 = x2 - 1.;
		cinv = 1. / hull(x2, 1.);
		r = 0.;
		xn = -1.;
		xn2 = -1.;
		for (i=1;  ; i++) {
			xn = -xn * x2m1; 
			xn2 = -xn2 * cinv * x2m1;
			tmp = xn2 / (T)i;
			if (mag(tmp) < std::numeric_limits<T>::epsilon()) {
				r += tmp;
				break;
			} else {
				r += xn / (T)i;
			}
		}

		r += constants<interval>::ln2() * p;

		if (round == -1) return r.lower();
		else return r.upper();
	}

	friend interval log(const interval& I) {
		return interval(log_point(I.lower(), -1), log_point(I.upper(), 1));
	}

	static interval log1p_origin(const T& x) {
		interval tmp;
		interval cinv;
		interval r;
		interval xn, xn2;
		int i;

		cinv = 1. / hull(x + interval(1.), 1.);
		r = 0.;
		xn = -1.;
		xn2 = -1.;
		for (i=1;  ; i++) {
			xn = -xn * x; 
			xn2 = -xn2 * cinv * x;
			tmp = xn2 / (T)(i);
			if (mag(tmp) < std::numeric_limits<T>::epsilon()) {
				r += xn2 / (T)(i);
				break;
			} else {
				r += xn / (T)(i);
			}
		}

		return r;
	}

	static T log1p_point(const T& x, int round) {
		interval tmp;

		if (x >= -(3. - 2. * std::sqrt(2.)) && x <= 3. - 2. * std::sqrt(2.)) {
			tmp = log1p_origin(x);
			if (round == -1) return tmp.lower();
			else return tmp.upper();
		} else {
			tmp = x + interval(1.);
			if (round == -1) {
				return log_point(tmp.lower(), -1);
			} else {
				return log_point(tmp.upper(), 1);
			}
		}
	}

	friend interval log1p(const interval& I) {
		return interval(log1p_point(I.lower(), -1), log1p_point(I.upper(), 1));
	}

	static interval sin_origin(const interval& I) {
		interval r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= I;
			y /= (T)i;
			if (mag(y) < std::numeric_limits<T>::epsilon()) {
				r += y * interval(-1., 1.);
				break;
			} else {
				if (i % 2 != 0) {
					if (i % 4 == 1) {
						r += y;
					} else {
						r -= y;
					}
				}
			}
		}

		return r;
	}

	static interval cos_origin(const interval& I) {
		interval r, y;
		int i;

		r = 1.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= I;
			y /= (T)i;
			if (mag(y) < std::numeric_limits<T>::epsilon()) {
				r += y * interval<T>(-1., 1.);
				break;
			} else {
				if (i % 2 == 0) {
					if (i % 4 == 0) {
						r += y;
					} else {
						r -= y;
					}
				}
			}
		}

		return r;
	}

	static interval sin_point(const interval& I) {
		const interval pi = constants<interval>::pi();
		const T mpi = mid(pi);

		if (I.lower() >= mpi) {
			return sin_point(I - pi * 2.);
		}

		if (I.upper() <= -mpi * 3. / 4.) {
			return -sin_origin(I + pi);
		}
		if (I.upper() <= -mpi * 0.5) {
			return -cos_origin(-pi * 0.5 - I);
		}
		if (I.upper() <= -mpi * 0.25) {
			return -cos_origin(I + pi * 0.5);
		}
		if (I.upper() <= 0.) {
			return -sin_origin(-I);
		}
		if (I.upper() <= mpi * 0.25) {
			return sin_origin(I);
		}
		if (I.upper() <= mpi * 0.5) {
			return cos_origin(pi * 0.5 - I);
		}
		if (I.upper() <= mpi * 3. / 4.) {
			return cos_origin(I - pi * 0.5);
		}
		return sin_origin(pi - I);
	}

	static interval cos_point(const interval& I) {
		const interval pi = constants<interval>::pi();
		const T mpi = mid(pi);

		if (I.lower() >= mpi) {
			return cos_point(I - pi * 2.);
		}

		if (I.upper() <= -mpi * 3. / 4.) {
			return -cos_origin(I + pi);
		}
		if (I.upper() <= -mpi * 0.5) {
			return -sin_origin(-pi * 0.5 - I);
		}
		if (I.upper() <= -mpi * 0.25) {
			return sin_origin(I + pi * 0.5);
		}
		if (I.upper() <= 0.) {
			return cos_origin(-I);
		}
		if (I.upper() <= mpi * 0.25) {
			return cos_origin(I);
		}
		if (I.upper() <= mpi * 0.5) {
			return sin_origin(pi * 0.5 - I);
		}
		if (I.upper() <= mpi * 3. / 4.) {
			return -sin_origin(I - pi * 0.5);
		}
		return -cos_origin(pi - I);
	}

	friend interval sin(const interval& I) {
		const interval pi = constants<interval>::pi();
		const interval pi2 = pi * 2.;

		T n;
		interval r, I2;

		using std::abs;
		if (abs(I.lower()) == std::numeric_limits<T>::infinity()) {
			return hull(-1. , 1.);
		}
		if (abs(I.upper()) == std::numeric_limits<T>::infinity()) {
			return hull(-1. , 1.);
		}

		I2 = I;
		while (I2.lower() <= -pi.upper() || I2.lower() >= pi.upper()) {
			using std::floor;
			n = floor((I2.lower() / pi2.lower()) + 0.5);
			I2 -= n * pi2;
		}

		if ((interval(I2.upper()) - I2.lower()).lower() >= pi2.upper()) {
			return interval(-1., 1.);
		}

		r = hull(sin_point(interval(I2.lower())), sin_point(interval(I2.upper())));

		if (subset(pi * 0.5, I2)) {
			r = hull(r, 1.);
		}

		if (subset(pi * 2.5, I2)) {
			r = hull(r, 1.);
		}

		if (subset(-pi * 0.5, I2)) {
			r = hull(r, -1.);
		}

		if (subset(pi * 1.5, I2)) {
			r = hull(r, -1.);
		}

		return intersect(r, interval(-1., 1.));
	}

	friend interval cos(const interval& I) {
		const interval pi = constants<interval>::pi();
		const interval pi2 = pi * 2.;

		T n;
		interval r, I2;

		using std::abs;
		if (abs(I.lower()) == std::numeric_limits<T>::infinity()) {
			return hull(-1. , 1.);
		}
		if (abs(I.upper()) == std::numeric_limits<T>::infinity()) {
			return hull(-1. , 1.);
		}

		I2 = I;
		while (I2.lower() <= -pi.upper() || I2.lower() >= pi.upper()) {
			using std::floor;
			n = floor((I2.lower() / pi2.lower()) + 0.5);
			I2 -= n * pi2;
		}

		if ((interval(I2.upper()) - I2.lower()).lower() >= pi2.upper()) {
			return interval(-1., 1.);
		}

		r = hull(cos_point(interval(I2.lower())), cos_point(interval(I2.upper())));

		if (in(0., I2)) {
			r = hull(r, 1.);
		}

		if (subset(pi2, I2)) {
			r = hull(r, 1.);
		}

		if (subset(-pi, I2)) {
			r = hull(r, -1.);
		}

		if (subset(pi, I2)) {
			r = hull(r, -1.);
		}

		if (subset(pi * 3., I2)) {
			r = hull(r, -1.);
		}

		return intersect(r, interval(-1., 1.));
	}

	static interval tan_point(const T& x) {
		return sin_point(interval(x)) / cos_point(interval(x));
	}

	friend interval tan(const interval& I) {
		const interval pi = constants<interval>::pi();
		const interval pih = pi * 0.5;

		T n;
		interval I2;

		using std::abs;
		if (abs(I.lower()) == std::numeric_limits<T>::infinity()) {
			return hull(-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
		}
		if (abs(I.upper()) == std::numeric_limits<T>::infinity()) {
			return hull(-std::numeric_limits<T>::infinity(), std::numeric_limits<T>::infinity());
		}

		I2 = I;
		while (I2.lower() <= -pih.upper() || I2.lower() >= pih.upper()) {
			using std::floor;
			n = floor((I2.lower() / pi.lower()) + 0.5);
			I2 -= n * pi;
		}

		if (I2.upper() >= pih.upper()) {
			return whole();
		}

		return interval(tan_point(I2.lower()).lower(), tan_point(I2.upper()).upper());
	}

	static interval atan_origin(const interval& I) {
		interval tmp;
		interval r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= I;
			tmp = y * interval(-1., 1.) / (T)i;
			if (mag(tmp) < std::numeric_limits<T>::epsilon()) {
				r += tmp;
				break;
			} else {
				if (i % 2 != 0) {
					if (i % 4 == 1) {
						r += y / (T)i;
					} else {
						r -= y / (T)i;
					}
				}
			}
		}

		return r;
	}

	static interval atan_point(const T& x) {
		const interval pi = constants<interval>::pi();

		interval I = interval(x);

		if (x < -(std::sqrt(2.) + 1.)) {
			return -pi * 0.5 - atan_origin(1. / I);
		}
		if (x < -(std::sqrt(2.) - 1.)) {
			return -pi * 0.25 + atan_origin((1. + I)/(1 - I));
		}
		if (x < (std::sqrt(2.) - 1.)) {
			return atan_origin(I);
		}
		if (x < (std::sqrt(2.) + 1.)) {
			return pi * 0.25 + atan_origin((I - 1.)/(I + 1.));
		}
		return pi * 0.5 - atan_origin(1. / I);
	}

	friend interval atan(const interval& I) {
		return interval(atan_point(I.lower()).lower(), atan_point(I.upper()).upper());
	}

	static interval asin_point(const T& x) {
		const interval pih = constants<interval>::pi() * 0.5;

		if (x == 1) return pih;
		if (x == -1) return -pih;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return atan(x / sqrt(1. - interval(x) * x));
		} else {
			if (x > 0.) {
				return atan(x / sqrt((1. + interval(x)) * (1. - x)));
			} else {
				return atan(x / sqrt((1. + x) * (1. - interval(x))));
			}
		}
	}

	friend interval asin(const interval<T>& I) {
		return interval(asin_point(I.lower()).lower(), asin_point(I.upper()).upper());
	}

	static interval pih_m_atan_point(const interval& I) {
		const interval pi = constants<interval>::pi();

		if (I.lower() < -(std::sqrt(2.) + 1.)) {
			return pi + atan_origin(1. / I);
		}
		if (I.lower() < -(std::sqrt(2.) - 1.)) {
			return pi * 0.75 - atan_origin((1. + I)/(1 - I));
		}
		if (I.lower() < (std::sqrt(2.) - 1.)) {
			return pi * 0.5 - atan_origin(I);
		}
		if (I.lower() < (std::sqrt(2.) + 1.)) {
			return pi * 0.25 - atan_origin((I - 1.)/(I + 1.));
		}
		return atan_origin(1. / I);
	}

	static interval acos_point(const T& x) {
		const interval pi = constants<interval>::pi();

		if (x == 1.) return interval(0.);
		if (x == -1.) return pi;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return pih_m_atan_point(x / sqrt(1. - interval(x) * x));
		} else {
			if (x > 0.) {
				return pih_m_atan_point(x / sqrt((1. + interval(x)) * (1. - x)));
			} else {
				return pih_m_atan_point(x / sqrt((1. + x) * (1. - interval(x))));
			}
		}
	}

	friend interval acos(const interval& I) {
		return interval(acos_point(I.upper()).lower(), acos_point(I.lower()).upper());
	}

	static interval atan2_point(const T& y, const T& x) {
		const interval pi = constants<interval>::pi();

		interval Ix = interval(x);
		interval Iy = interval(y);

		if (y <= x && y > -x) {
			return atan(Iy / Ix);
		}
		if (y > x && y > -x) {
			return pi * 0.5 - atan(Ix / Iy);
		}
		if (y > x && y <= -x) {
			if (y >= 0.) {
				return pi + atan(Iy / Ix);
			} else {
				return -pi + atan(Iy / Ix);
			}
		}
		return -pi * 0.5 - atan(Ix / Iy);
	}

	friend interval atan2(const interval& Iy, const interval& Ix) {
		const interval pi = constants<interval>::pi();

		if (zero_in(Ix)) {
			if (zero_in(Iy)) {
				return interval(-pi.upper(), pi.upper());
			} else {
				if (Iy > 0.) {
					return interval(atan2_point(Iy.lower(), Ix.upper()).lower(), atan2_point(Iy.lower(), Ix.lower()).upper());
				} else {
					return interval(atan2_point(Iy.upper(), Ix.lower()).lower(), atan2_point(Iy.upper(), Ix.upper()).upper());
				}
			}
		} else {
			if (zero_in(Iy)) {
				if (Ix > 0.) {
					return interval(atan2_point(Iy.lower(), Ix.lower()).lower(), atan2_point(Iy.upper(), Ix.lower()).upper());
				} else {
					if (Iy.lower() < 0)
						return interval(atan2_point(Iy.upper(), Ix.upper()).lower(), (pi * 2. + atan2_point(Iy.lower(), Ix.upper())).upper());
					else {
						return interval(atan2_point(Iy.upper(), Ix.upper()).lower(), (atan2_point(Iy.lower(), Ix.upper())).upper());
					}
				}
			} else {
				if (Ix > 0.) {
					if (Iy > 0.) {
						return interval(atan2_point(Iy.lower(), Ix.upper()).lower(), atan2_point(Iy.upper(), Ix.lower()).upper());
					} else {
						return interval(atan2_point(Iy.lower(), Ix.lower()).lower(), atan2_point(Iy.upper(), Ix.upper()).upper());
					}
				} else {
					if (Iy > 0.) {
						return interval(atan2_point(Iy.upper(), Ix.upper()).lower(), atan2_point(Iy.lower(), Ix.lower()).upper());
					} else {
						return interval(atan2_point(Iy.upper(), Ix.lower()).lower(), atan2_point(Iy.lower(), Ix.upper()).upper());
					}
				}
			}
		}
	}

	static interval sinh_origin(const T& x) {
		// cosh(0.5) = (exp(0.5)+exp(-0.5))/2
		const interval exph = sqrt(constants<interval>::e());
		const interval coshh = (exph + 1. / exph) * 0.5;

		interval tmp;
		interval r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= x;
			y /= (T)i;
			tmp = y * interval(-coshh.upper(), coshh.upper());
			if (mag(tmp) < std::numeric_limits<T>::epsilon()) {
				r += tmp;
				break;
			} else {
				if (i % 2 != 0) {
					r += y;
				}
			}
		}

		return r;
	}

	static interval sinh_point(const T& x) {
		if (x >= -0.5 && x <= 0.5) {
			return sinh_origin(x);
		} else {
			interval tmp;
			tmp = exp_point(x);
			return (tmp - 1./tmp) * 0.5;
		}
	}

	friend interval sinh(const interval& I) {
		return interval(sinh_point(I.lower()).lower(), sinh_point(I.upper()).upper());
	}

	static interval cosh_point(const T& x) {
		interval tmp;
		tmp = exp_point(x);
		return (tmp + 1./tmp) * 0.5;
	}

	friend interval cosh(const interval& I) {
		interval r;
		r = hull(cosh_point(I.lower()), cosh_point(I.upper()));
		if (zero_in(I)) {
			r = hull(r, 1.);
		}
		return r;
	}

	static interval tanh_point(const T& x) {
		return sinh_point(x) / cosh_point(x);
	}

	friend interval tanh(const interval& I) {
		return interval(tanh_point(I.lower()).lower(), tanh_point(I.upper()).upper());
	}

	static interval asinh_point(const T& x) {
		if (x < -0.5) {
			return -log(-x + sqrt(1. + interval(x) * x));
		} else if (x <= 0.5) {
			return log1p((1. + x / (1. + sqrt(1. + interval(x) * x))) * x);
		} else {
			return log(x + sqrt(1. + interval(x) * x));
		}
	}

	friend interval asinh(const interval& I) {
		return interval(asinh_point(I.lower()).lower(), asinh_point(I.upper()).upper());
	}

	static interval acosh_point(const T& x) {
		if (x == 1.) {
			return interval(0.);
		} else if (x <= 1.5) {
			interval y(x - 1.);
			return log1p(y + sqrt(y * (interval(x) + 1.)));
		} else {
			return log(x + sqrt(interval(x) * x - 1.));
		}
	}

	friend interval acosh(const interval& I) {
		return interval(acosh_point(I.lower()).lower(), acosh_point(I.upper()).upper());
	}

	static interval atanh_point(const T& x) {
		if (x < -0.5) {
			return 0.5 * log((1. + x) / (1. - interval(x)));
		} else if (x <= 0.5) {
			return 0.5 * log1p(2. * x / (1. - interval(x)));
		} else {
			return 0.5 * log((1. + interval(x)) / (1. - x));
		}
	}

	friend interval atanh(const interval& I) {
		return interval(atanh_point(I.lower()).lower(), atanh_point(I.upper()).upper());
	}
};


template <class T> struct constants< interval<T> > {
	static interval<T> pi() {
		// static is used so that string is evaluated only "one time"
		static const interval<T> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<T> e() {
		static const interval<T> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<T> ln2() {
		static const interval<T> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}
	static interval<T> str(const std::string& s) {
		return interval<T>(s, s);
	}
	static interval<T> str(const std::string& s1, const std::string& s2) {
		return interval<T>(s1, s2);
	}
};

} // namespace kv

#endif // INTERVAL_HPP
