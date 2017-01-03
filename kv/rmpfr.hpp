/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RMPFR_HPP
#define RMPFR_HPP


namespace kv {

template <int N> struct rop < mpfr<N> > {

	static mpfr<N> add_up(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_add(z.a, x.a, y.a, MPFR_RNDU);
		return z;
	}

	static mpfr<N> add_down(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_add(z.a, x.a, y.a, MPFR_RNDD);
		return z;
	}

	static mpfr<N> sub_up(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_sub(z.a, x.a, y.a, MPFR_RNDU);
		return z;
	}

	static mpfr<N> sub_down(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_sub(z.a, x.a, y.a, MPFR_RNDD);
		return z;
	}

	static mpfr<N> mul_up(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_mul(z.a, x.a, y.a, MPFR_RNDU);
		return z;
	}

	static mpfr<N> mul_down(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_mul(z.a, x.a, y.a, MPFR_RNDD);
		return z;
	}

	static mpfr<N> div_up(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_div(z.a, x.a, y.a, MPFR_RNDU);
		return z;
	}

	static mpfr<N> div_down(const mpfr<N>& x, const mpfr<N>& y) {
		mpfr<N> z;
		mpfr_div(z.a, x.a, y.a, MPFR_RNDD);
		return z;
	}

	static mpfr<N> sqrt_up(const mpfr<N>& x){
		mpfr<N> z;
		mpfr_sqrt(z.a, x.a, MPFR_RNDU);
		return z;
	}

	static mpfr<N> sqrt_down(const mpfr<N>& x) {
		mpfr<N> z;
		mpfr_sqrt(z.a, x.a, MPFR_RNDD);
		return z;
	}

	static void begin() {
	}

	static void end() {
	}

	static void print_up(const mpfr<N>& x, std::ostream& s) {
		char format;
		char fmt[100];
		char *str;

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

		sprintf(fmt, "%%.%dRU%c", (int)s.precision(), format);
		mpfr_asprintf(&str, fmt, x.a);
		s << str;
		mpfr_free_str(str);
	}

	static void print_down(const mpfr<N>& x, std::ostream& s) {
		char format;
		char fmt[100];
		char *str;

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

		sprintf(fmt, "%%.%dRD%c", (int)s.precision(), format);
		mpfr_asprintf(&str, fmt, x.a);
		s << str;
		mpfr_free_str(str);
	}

	static mpfr<N> fromstring_up(const std::string& s) {
		mpfr<N> r;
		mpfr_set_str(r.a, s.c_str(), 10, MPFR_RNDU);
		return r;
	}

	static mpfr<N> fromstring_down(const std::string& s) {
		mpfr<N> r;
		mpfr_set_str(r.a, s.c_str(), 10, MPFR_RNDD);
		return r;
	}
};

template <int N> struct constants< interval< mpfr<N> > > {
	static interval< mpfr<N> > pi() {
		static interval< mpfr<N> > tmp(0);
		if (tmp.lower() != 0) return tmp;
		mpfr_const_pi(tmp.lower().a, MPFR_RNDD);
		mpfr_const_pi(tmp.upper().a, MPFR_RNDU);
		return tmp;
	}

	static interval< mpfr<N> > e() {
		static interval< mpfr<N> > tmp(0);
		if (tmp.lower() != 0) return tmp;
		mpfr<N> one(1);
		mpfr_exp(tmp.lower().a, one.a, MPFR_RNDD);
		mpfr_exp(tmp.upper().a, one.a, MPFR_RNDU);
		return tmp;
	}

	static interval< mpfr<N> > ln2() {
		static interval< mpfr<N> > tmp(0);
		if (tmp.lower() != 0) return tmp;
		mpfr_const_log2(tmp.lower().a, MPFR_RNDD);
		mpfr_const_log2(tmp.upper().a, MPFR_RNDU);
		return tmp;
	}

	static interval< mpfr<N> > str(const std::string& s) {
		return interval< mpfr<N> >(s, s);
	}

	static interval< mpfr<N> > str(const std::string& s1, const std::string& s2) {
		return interval< mpfr<N> >(s1, s2);
	}
};

} // namespace kv

#endif // RMPFR_HPP
