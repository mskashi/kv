#include <cstdio>
#include <iostream>
#include <string>
#include <cctype>
#include <vector>
#include <list>
#include <limits>
#include <cmath>
#include <cstdlib>
#include "dd.hpp"

inline int get_sign_double(double x)
{
	if (x == 0.) {
		x = 1. / x;
	}

	if (x > 0.) return 1;
	else return -1;
}

inline int get_exponent(double x)
{
	int i;

	if (x >= ldexp(1., 1023)) return 1023;
	if (x < ldexp(1., -1074)) return -1075;

	frexp(x, &i);

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

inline std::string ddtostring(kv::dd x, int precision = 34, char format = 'g', int mode = 0)
{
	int i, j;
	int sign, sign2, ex1, ex2;
	double absx1, absx2;

	sign = get_sign_double(x.a1);
	absx1 = fabs(x.a1);

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

	// get x.a1 to buf

	ex1 = get_exponent(absx1);

	// add 1-byte margin to add buf2
	bool buf[1023 - (-1074) + 2];
	int offset = 1074;
	int emax, emin;
	double dtmp, dtmp2;

	dtmp = absx1;
	dtmp2 = ldexp(1., ex1);

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

	// get x.a2 to buf2 and add it to buf

	bool buf2[1023 - (-1074) + 1];
	int emax2, emin2, s;
	int carry, tmp;

	sign2 = get_sign_double(x.a2);
	absx2 = fabs(x.a2);

	if (absx2 != 0.) {
		ex2 = get_exponent(absx2);
		dtmp = absx2;
		dtmp2 = ldexp(1., ex2);

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
			carry = (int)floor(tmp / 2.);
			buf[offset + i] = tmp % 2;
		}
		for (i=emax2+1; emax; i++) {
			if (carry == 0) break;
			// NOTICE: tmp may become negative
			tmp = buf[offset + i] + carry;
			carry = (int)floor(tmp / 2.);
			buf[offset + i] = tmp % 2;
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

		while (buf[offset + emax] == 0 and emax >= 0) {
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

		while (buf[offset + emin] == 0 && emin < 0) {
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

	std::string result_str;
	char stmp[100];

	if (sign == 1) {
		result_str = "";
	} else {
		result_str = "-";
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
					if (result[i] != 10) break;
					result[offset2 + i] = 0;
				}
				if (result[offset2 + result_max] == 0) {
					result_max--;
				}
			}
		}

		// delete zeros of tail
		while (result[offset2 + result_min] == 0 && result_min < 0) {
			result_min++;
		}

		// make result string
		for (i=result_max; i>=result_min; i--) {
			if (i == -1) result_str += ".";
			sprintf(stmp, "%d", result[offset2 + i]);
			result_str += stmp;
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
					if (result[i] != 10) break;
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
			if (i == result_max -1) result_str += ".";
			sprintf(stmp, "%d", result[offset2 + i]);
			result_str += stmp;
		}
		sprintf(stmp, "e%+03d", result_max);
		result_str += stmp;

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
					if (result[i] != 10) break;
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
			while (result[offset2 + result_min] == 0 && result_min < 0) {
				result_min++;
			}

			if (result_max < 0) {
				result_max = 0;
			}

			// make result string
			for (i=result_max; i>=result_min; i--) {
				if (i == -1) result_str += ".";
				sprintf(stmp, "%d", result[offset2 + i]);
				result_str += stmp;
			}

		} else {
			// use 'e' like format

			// delete zeros of tail
			while (result[offset2 + result_min] == 0) {
				result_min++;
			}

			// make result string
			for (i=result_max; i>=result_min; i--) {
				if (i == result_max -1) result_str += ".";
				sprintf(stmp, "%d", result[offset2 + i]);
				result_str += stmp;
			}
			sprintf(stmp, "e%+03d", result_max);
			result_str += stmp;
		}

	} else if (format == 'a') {
		// make result string
		for (i=result_max; i>=result_min; i--) {
			if (i == -1) result_str += ".";
			sprintf(stmp, "%d", result[offset2 + i]);
			result_str += stmp;
		}
	}

	return result_str;
}


int main()
{
	kv::dd x;
	int i;
	std::cout.precision(20);

	x = kv::dd(1.) / 3.;

	for (i=0; i<10; i++) x /= 10.;

	for (i=0; i<20; i++) {
		std::cout << "e: " << ddtostring(x, 34, 'e', 0) << "\n";
		std::cout << "f: " << ddtostring(x, 34, 'f', 0) << "\n";
		std::cout << "g: " << ddtostring(x, 34, 'g', 0) << "\n";
		std::cout << "a: " << ddtostring(x, 34, 'a', 0) << "\n";
		x *= 10.;
	}
}
