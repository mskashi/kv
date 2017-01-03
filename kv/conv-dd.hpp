/*
 * Copyright (c) 2013-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef CONV_DD_HPP
#define CONV_DD_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cctype>
#include <vector>
#include <list>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <algorithm>


namespace kv {

struct conv_dd {

	static int get_sign_double(double x) {
		if (x == 0.) {
			x = 1. / x;
		}

		if (x > 0.) return 1;
		else return -1;
	}

	static int get_exponent(double x) {
		int i;

		if (x >= std::ldexp(1., 1023)) return 1023;
		if (x < std::ldexp(1., -1074)) return -1075;

		std::frexp(x, &i);

		return i - 1;
	}

	// convert double-double number to string
	// mode == -1 : down
	// mode ==  0 : nearest
	// mode ==  1 : up
	// format == 'e' : like %e of printf
	// format == 'f' : like %f of printf
	// format == 'g' : like %g of printf
	// format == 'a' : print all digits with no rounding

	static std::string ddtostring(double x1, double x2, int precision = 34, char format = 'g', int mode = 0) {
		int i, j;
		int sign, sign2, ex1, ex2;
		double absx1, absx2;

		if (x1 != x1 || x2 != x2) return "nan";

		sign = get_sign_double(x1);
		absx1 = std::fabs(x1);

		if (absx1 == 0.) {
			if (sign == -1) {
				return "-0";
			} else {
				return "0";
			}
		}

		if (absx1 == std::numeric_limits<double>::infinity()) {
			if (sign == -1) {
				return "-inf";
			} else {
				return "inf";
			}
		}

		// get x1 to buf

		ex1 = get_exponent(absx1);

		// add 1-byte margin to add buf2
		bool buf[1023 - (-1074) + 2];
		int offset = 1074;
		int emax, emin;
		double dtmp, dtmp2;

		dtmp = absx1;
		dtmp2 = std::ldexp(1., ex1);

		for (i=0; i<=52; i++) {
			if (dtmp >= dtmp2) {
				buf[offset + ex1 - i] = 1;
				dtmp -= dtmp2;
			} else {
				buf[offset + ex1 - i] = 0;
			}
			if (dtmp == 0) {
				emax = ex1;
				emin = ex1 - i;
				break;
			}
			dtmp2 /= 2.;
		}

		// get x2 to buf2 and add it to buf

		bool buf2[1023 - (-1074) + 1];
		int emax2, emin2, s;
		int carry, tmp;

		sign2 = get_sign_double(x2);
		absx2 = std::fabs(x2);

		if (absx2 != 0.) {
			ex2 = get_exponent(absx2);
			dtmp = absx2;
			dtmp2 = std::ldexp(1., ex2);

			for (i=0; i<=52; i++) {
				if (dtmp >= dtmp2) {
					buf2[offset + ex2 - i] = 1;
					dtmp -= dtmp2;
				} else {
					buf2[offset + ex2 - i] = 0;
				}
				if (dtmp == 0) {
					emax2 = ex2;
					emin2 = ex2 - i;
					break;
				}
				dtmp2 /= 2.;
			}

			if (sign == sign2)  s = 1;
			else s = -1;

			if (emin > emin2) {
				for (i=emin2; i<=emin-1; i++) {
					buf[offset + i] = 0;
				}
				emin = emin2;
			}
			emax++;
			buf[offset + emax] = 0;

			carry = 0;
			for (i=emin2; i<=emax2; i++) {
				// NOTICE: tmp may become negative
				tmp = buf[offset + i] + s * buf2[offset + i] + carry;
				carry = (int)std::floor(tmp / 2.);
				buf[offset + i] = tmp - carry * 2;
			}
			for (i=emax2+1; i<=emax; i++) {
				if (carry == 0) break;
				// NOTICE: tmp may become negative
				tmp = buf[offset + i] + carry;
				carry = (int)std::floor(tmp / 2.);
				buf[offset + i] = tmp - carry * 2;
			}
			while (buf[offset + emax] == 0) {
				emax--;
			}
		}


		if (emin > 0) {
			for (i=0; i<emin; i++) {
				buf[offset + i] = 0;
			}
			emin = 0;
		}

		if (emax < 0) {
			for (i=emax + 1; i<=0; i++) {
				buf[offset + i] = 0;
			}
			emax = 0;
		}

		std::list<int> result1, result2;
		int result_max, result_min, m, pm;

		result_max = -1;

		while (emax >= 0) {
			if (emax >= 17) m= 5;
			else if (emax >= 14) m = 4;
			else if (emax >= 10) m = 3;
			else if (emax >= 7) m = 2;
			else  m = 1;

			pm = 1;
			for (i=0; i<m; i++) pm *= 10;

			carry = 0;
			for (i=emax; i>=0; i--) {
				tmp = carry * 2 + buf[offset + i];
				buf[offset + i] = tmp / pm;
				carry = tmp % pm;
			}

			for (i=0; i<m; i++) {
				result_max++;
				result1.push_back(carry  % 10);
				carry /= 10;
			}

			while (emax >= 0 && buf[offset + emax] == 0) {
				emax--;
			}
		}

		result_min = 0;

		while (emin < 0) {
			m = std::min(8, -emin);
			pm = 1;
			for (i=0; i<m; i++) pm *= 10;

			carry = 0;
			for (i=emin; i<=-1; i++) {
				tmp = buf[offset + i] * pm + carry;
				buf[offset + i] = tmp % 2;
				carry = tmp / 2;
			}

			for (i=0; i<m; i++) {
				result_min--;
				pm /= 10;
				result2.push_back(carry / pm);
				carry %= pm;
			}

			while (emin < 0 && buf[offset + emin] == 0) {
				emin++;
			}
		}

		std::vector<int> result;
		int offset2;

		// add 1byte margin to both ends of array
		result.resize(result_max - result_min + 1 + 2);
		offset2 = - result_min + 1;
		for (i=0; i<=result_max; i++) {
			result[offset2 + i] = result1.front();
			result1.pop_front();
		}
		for (i=-1; i>=result_min; i--) {
			result[offset2 + i] = result2.front();
			result2.pop_front();
		}

		#if 0
		for (i=result_min; i<=result_max; i++) {
			std::cout << i << ':' << result[offset2 + i] << "\n";
		}
		#endif

		std::ostringstream result_str;

		if (sign == 1) {
			result_str << "";
		} else {
			result_str << "-";
		}

		if (format == 'f') {
			// round to precision after decimal point
			if (-(precision+1) >= result_min) {
				result_min = -precision;
				tmp = result[offset2 + result_min - 1];
				if ((mode == 1 && sign == 1) || (mode == -1 && sign == -1) || (mode == 0 && tmp >= 5)) {
					result[offset2 + result_max + 1] = 0;
					result_max++;
					for (i=result_min; i<=result_max; i++) {
						result[offset2 + i]++;
						if (result[offset2 + i] != 10) break;
						result[offset2 + i] = 0;
					}
					if (result[offset2 + result_max] == 0) {
						result_max--;
					}
				}
			}

			// delete zeros of tail
			while (result_min < 0 && result[offset2 + result_min] == 0) {
				result_min++;
			}

			// make result string
			for (i=result_max; i>=result_min; i--) {
				if (i == -1) result_str << ".";
				result_str << result[offset2 + i];
			}

		} else if (format == 'e') {
			// delete zeros of head
			while (result[offset2 + result_max] == 0) {
				result_max--;
			}

			// round to precision
			if (result_max-precision-1 >= result_min) {
				result_min = result_max - precision;
				tmp = result[offset2 + result_min - 1];
				if ((mode == 1 && sign == 1) || (mode == -1 && sign == -1) || (mode == 0 && tmp >= 5)) {
					result[offset2 + result_max + 1] = 0;
					result_max++;
					for (i=result_min; i<=result_max; i++) {
						result[offset2 + i]++;
						if (result[offset2 + i] != 10) break;
						result[offset2 + i] = 0;
					}
					if (result[offset2 + result_max] == 0) {
						result_max--;
					} else {
						result_min++;
					}
				}
			}

			// delete zeros of tail
			while (result[offset2 + result_min] == 0) {
				result_min++;
			}

			// make result string
			for (i=result_max; i>=result_min; i--) {
				if (i == result_max -1) result_str << ".";
				result_str << result[offset2 + i];
			}
			// sprintf(stmp, "e%+03d", result_max);
			result_str << "e" << std::internal << std::showpos << std::setfill('0') << std::setw(3) << result_max << std::noshowpos;

		} else if (format == 'g') {
			// delete zeros of head
			while (result[offset2 + result_max] == 0) {
				result_max--;
			}

			// round to precision
			if (result_max-precision >= result_min) {
				result_min = result_max - precision + 1;
				tmp = result[offset2 + result_min - 1];
				if ((mode == 1 && sign == 1) || (mode == -1 && sign == -1) || (mode == 0 && tmp >= 5)) {
					result[offset2 + result_max + 1] = 0;
					result_max++;
					for (i=result_min; i<=result_max; i++) {
						result[offset2 + i]++;
						if (result[offset2 + i] != 10) break;
						result[offset2 + i] = 0;
					}
					if (result[offset2 + result_max] == 0) {
						result_max--;
					} else {
						result_min++;
					}
				}
			}

			if (-4 <= result_max && result_max <= precision -1) {
				// use 'f' like format

				// delete zeros of tail
				while (result_min < 0 && result[offset2 + result_min] == 0) {
					result_min++;
				}

				if (result_max < 0) {
					result_max = 0;
				}

				// make result string
				for (i=result_max; i>=result_min; i--) {
					if (i == -1) result_str << ".";
					result_str << result[offset2 + i];
				}

			} else {
				// use 'e' like format

				// delete zeros of tail
				while (result[offset2 + result_min] == 0) {
					result_min++;
				}

				// make result string
				for (i=result_max; i>=result_min; i--) {
					if (i == result_max -1) result_str << ".";
					result_str << result[offset2 + i];
				}
				// sprintf(stmp, "e%+03d", result_max);
				result_str << "e" << std::internal << std::showpos << std::setfill('0') << std::setw(3) << result_max << std::noshowpos;
			}

		} else if (format == 'a') {
			// make result string
			for (i=result_max; i>=result_min; i--) {
				if (i == -1) result_str << ".";
				result_str << result[offset2 + i];
			}
		}

		return result_str.str();
	}


	static void ignore_space(std::string& s) {
		int p = 0;
		while (p < s.size() && isspace(s[p])) p++;
		s = s.substr(p);
	}

	static int get_sign(std::string& s) {
		int r, p;

		p = 0;
		r = 1;
		if (p < s.size() && s[p] == '-') {
			r = -1;
			p++;
		} else if (p < s.size() && s[p] == '+') {
			r = 1;
			p++;
		}

		s = s.substr(p);
		return r;
	}

	static std::string get_number(std::string& s) {
		std::string r;
		int p;

		p = 0;
		r = "";
		while (p < s.size() && isdigit(s[p])) {
			r += s[p];
			p++;
		}

		s = s.substr(p);
		return r;
	}

	// convert string to double-double number
	// mode == -1 : down
	// mode ==  0 : nearest
	// mode ==  1 : up

	// REMARK: if fast is true, we convert only leading 107bit. 
	//  Consequently, the second part of double-double number 
	//  may not achieve "best" precision.

	static void stringtodd(std::string s, double& x1, double& x2, int mode = 0, bool fast = false) {
		int i, j, tmp;
		bool flag;
		int sign, e10, esign;
		std::string num1_s, num2_s, nume_s;

		ignore_space(s);
		sign = get_sign(s);
		num1_s = get_number(s);

		if (0 < s.size() && s[0] == '.') {
			s = s.substr(1);
			num2_s = get_number(s);
		}

		if (0 < s.size() && (s[0] == 'e' || s[0] == 'E')) {
			s = s.substr(1);
			esign = get_sign(s);
			nume_s = get_number(s);
			e10 = esign * atoi(nume_s.c_str());
		} else {
			e10 = 0;
		}

		// delete 0s from the head of num1_s
		while (0 < num1_s.size() && num1_s[0] == '0') {
			num1_s = num1_s.substr(1);
		}

		// delete 0s from the tail of num2_s
		while (0 < num2_s.size() && num2_s[num2_s.size() - 1] == '0') {
			num2_s = num2_s.substr(0, num2_s.size() - 1);
		}

		// set table and offset
		// |x| = \sum_{table_min}^{table_max} table[offset + i] * 10^i
		int table_max, table_min, offset;
		std::vector<int> table;

		table_max = num1_s.size() - 1 + e10;
		table_min = e10 - num2_s.size();
		table.resize(table_max - table_min + 1);
		offset = - table_min;

		for (i=0; i<num1_s.size(); i++) {
			table[offset + num1_s.size() - 1 - i + e10] = num1_s[i] - '0';
		}

		for (i=0; i<num2_s.size(); i++) {
			table[offset - i - 1 + e10] = num2_s[i] - '0';
		}

		// extend table
		if (table_min > 0) {
			tmp = table.size();
			table.resize(tmp + table_min);
			for (i=tmp-1; i>=0; i--) {
				table[i + table_min] = table[i];
			}
			for (i=0; i<table_min; i++) {
				table[i] = 0;
			}
			offset += table_min;
			table_min = 0;
		}

		if (table_max < -1) {
			tmp = table.size();
			table.resize(tmp + (-1-table_max));
			for (i=0; i<(-1-table_max); i++) {
				table[tmp + i] = 0;
			}
			table_max = -1;
		}

		#if 0
		for (i=table_max; i>=table_min; i--) {
			std::cout << i << ':' << table[offset + i] << "\n";
		}
		#endif

		// convert decimal number to binary number
		// set result and offset2
		// |x| = \sum_{result_min}^{reuslt_max} result[offset2 + i] * 2^i

		int result_min, result_max, m, pm, carry, carry2;
		std::list<bool> result1, result2;

		// integer part

		result_max = -1;

		while (table_max >= 0) {
			if (table_max >= 5) m = 16;
			else if (table_max >= 4) m = 13;
			else if (table_max >= 3) m = 9;
			else if (table_max >= 2) m = 6;
			else if (table_max >= 1) m = 3;
			else m = 1;
			pm = 1 << m;
			
			carry = 0;
			for (i=table_max; i>=0; i--) {
				tmp = carry * 10 + table[offset + i];
				table[offset + i] = tmp / pm;
				carry = tmp % pm;
			}
			for (i=0; i<m; i++) {
				result_max++;
				result1.push_back(carry % 2);
				carry = carry / 2;
			}
			while (table_max >= 0 && table[offset + table_max] == 0) {
				table_max--;
			}
		}

		// fraction part

		//  flag means whether most significant bit already found or not
		if (result_max >= 0) flag = true;
		else flag = false;

		result_min = 0;

		while (table_min < 0) {
			if (fast) {
				tmp = 106 - (result_max - result_min);
			} else {
				tmp = result_min + 1075;
			}
			if (flag && tmp <= 0) break;
			if (!flag) {
				m = 16;
			} else {
				m = std::min(16, tmp);
			}
			pm = 1 << m;

			carry = 0;
			for (i=table_min; i<=-1; i++) {
				tmp = table[offset + i] * pm + carry;
				table[offset + i] = tmp % 10;
				carry = tmp / 10;
			}

			for (i=0; i<m; i++) {
				result_min--;
				pm /= 2;
				carry2 = carry / pm;
				carry = carry % pm;

				if (flag) {
					result2.push_back(carry2);
				} else {
					if (carry2 != 0) {
						result2.push_back(carry2);
						result_max = result_min;
						flag = true;
					}
				}
			}

			while (table_min < 0 && table[offset + table_min] == 0) {
				table_min++;
			}
		}

		// append integer and fraction part

		std::vector<bool> result;
		int offset2;

		result.resize(result_max - result_min + 1);
		offset2 = - result_min;
		for (i=0; i<=result_max; i++) {
			result[offset2 + i] = result1.front();
			result1.pop_front();
		}
		for (i=std::min(-1, result_max); i>=result_min; i--) {
			result[offset2 + i] = result2.front();
			result2.pop_front();
		}

		#if 0
		for (i=result_max; i>=result_min; i--) {
			printf("%d %d\n", i, result[offset2 + i]);
		}
		#endif

		// convert binary to double double number

		double dtmp;

		if (result_max > 1023) {
			if ((sign == 1 && mode == -1) || (sign == -1 && mode == 1)) {
				x1 = sign * (std::numeric_limits<double>::max)();
				x2 = std::ldexp(x1, -54);
				return;
			}
			dtmp = sign * std::numeric_limits<double>::infinity();
			x1 = dtmp;
			x2 = dtmp;
			return;
		}

		if (result_max < -1075) {
			if ((sign == 1 && mode == 1) || (sign == -1 && mode == -1)) {
				x1 = sign * std::numeric_limits<double>::denorm_min();
				x2 = sign * 0.;
				return;
			}
			dtmp = sign * 0.;
			x1 = dtmp;
			x2 = dtmp;
			return;
		}

		double r, r2;
		int result_max2;
		int msb;

		r = 0.;
		flag = false; // roundup first part or not
		result_max2 = result_min - 1;
		for (i=result_max; i >= result_min; i--) {
			if (result_max - i == 53 || i == -1075) {
				if (sign == 1) {
					if (result[offset2 + i] == 0) {
					} else {
						r += std::ldexp(1., i+1);
						flag = true;
					}
				} else {
					if (result[offset2 + i] == 0) {
					} else {
						r += std::ldexp(1., i+1);
						flag = true;
					}
				}
				result_max2 = i;
				break;
			}
			r += std::ldexp((double)result[offset2 + i], i);
		}

		if (flag) {
			r2 = - std::ldexp(1., result_max2 + 1);
		} else {
			r2 = 0.;
		}

		msb = result_min - 1; // outside result bits

		for (i=result_max2; i >= result_min; i--) {
			if (fast) {
				tmp = result_max2 - i;
			} else {
				tmp = msb - i;
			}
			if (tmp == 53 || i == -1075) {
				if (sign == 1) {
					if (mode == -1) {
					} else if (mode == 0) {
						if (result[offset2 + i] == 0) {
						} else {
							r2 += std::ldexp(1., i+1);
						}
					} else {
						r2 += std::ldexp(1., i+1);
					}
				} else {
					if (mode == -1) {
						r2 += std::ldexp(1., i+1);
					} else if (mode == 0) {
						if (result[offset2 + i] == 0) {
						} else {
							r2 += std::ldexp(1., i+1);
						}
					} else {
					}
				}
				break;
			}
			tmp = result[offset2 + i];
			r2 += std::ldexp((double)tmp, i);
			if (msb == result_min - 1 && ((flag && tmp == 0) || (!flag && tmp == 1))) {
				msb = i;
			}
		}

		x1 = sign * r;
		x2 = sign * r2;
		return;
	}
};

} // namespace kv


#endif // CONV_DD_HPP
