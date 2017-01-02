#include <iostream>
#include <string>
#include <cctype>
#include <vector>
#include <list>
#include <limits>
#include "dd.hpp"

inline void ignore_space(std::string& s)
{
	int p = 0;
	while (isspace(s[p])) p++;
	s = s.substr(p);
}

inline int get_sign(std::string& s)
{
	int r, p;

	p = 0;
	r = 1;
	if (s[p] == '-') {
		r = -1;
		p++;
	} else if (s[p] == '+') {
		r = 1;
		p++;
	}

	s = s.substr(p);
	return r;
}

inline std::string get_number(std::string& s)
{
	std::string r;
	int p;

	p = 0;
	r = "";
	while (isdigit(s[p])) {
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

inline kv::dd stringtodd(std::string s, int mode = 0, bool fast = false)
{
	int i, j, tmp;
	bool flag;
	int sign, e10, esign;
	std::string num1_s, num2_s, nume_s;

	ignore_space(s);
	sign = get_sign(s);
	num1_s = get_number(s);

	if (s[0] == '.') {
		s = s.substr(1);
		num2_s = get_number(s);
	}

	if (s[0] == 'e' || s[0] == 'E') {
		s = s.substr(1);
		esign = get_sign(s);
		nume_s = get_number(s);
		e10 = esign * atoi(nume_s.c_str());
	} else {
		e10 = 0;
	}

	// delete 0s from the head of num1_s
	while (num1_s[0] == '0') {
		num1_s = num1_s.substr(1);
	}

	// delete 0s from the tail of num2_s
	while (num2_s[num2_s.size() - 1] == '0') {
		num2_s = num2_s.substr(0, num2_s.size() - 1);
	}

	// set table and offset
	// |x| = \sum_{table_min}^{table_max} table[offset + i] * 10^i
	int table_max, table_min, offset;
	std::vector<int> table;

	table_max = num1_s.size() - 1 + e10;
	table_min = - num2_s.size() + e10;
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
		while (table[offset + table_max] == 0 && table_max >= 0) {
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

		while (table[offset + table_min] == 0 && table_min < 0) {
			table_min--;
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
		dtmp = sign * std::numeric_limits<double>::infinity();
		return kv::dd(dtmp, dtmp);
	}

	if (result_max < -1075) {
		dtmp = sign * 0.;
		return kv::dd(dtmp, dtmp);
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
					r += ldexp(1., i+1);
					flag = 1;
				}
			} else {
				if (result[offset2 + i] == 0) {
				} else {
					r += ldexp(1., i+1);
					flag = 1;
				}
			}
			result_max2 = i;
			break;
		}
		r += ldexp((double)result[offset2 + i], i);
	}

	if (flag) {
		r2 = - ldexp(1., result_max2 + 1);
	} else {
		r2 = 0.;
	}

	msb = std::numeric_limits<int>::min();

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
						r2 += ldexp(1., i+1);
					}
				} else {
					r2 += ldexp(1., i+1);
				}
			} else {
				if (mode == -1) {
					r2 += ldexp(1., i+1);
				} else if (mode == 0) {
					if (result[offset2 + i] == 0) {
					} else {
						r2 += ldexp(1., i+1);
					}
				} else {
				}
			}
			break;
		}
		tmp = result[offset2 + i];
		r2 += ldexp((double)tmp, i);
		if (msb == std::numeric_limits<int>::min() && ((flag == true && tmp == 0) || (flag == false && tmp == 1))) {
			msb = i;
		}
	}

	return kv::dd(sign * r, sign * r2);
}


int main()
{
	std::cout << std::cout.precision() << "\n";
	std::cout.precision(200);
	std::cout << std::cout.precision() << "\n";
	std::cout << (std::cout.flags() & std::cout.scientific) << "\n";
	std::cout << std::cout.scientific << "\n";
	std::cout.setf(std::cout.scientific);
	std::cout << (std::cout.flags() & std::cout.scientific) << "\n";
	std::cout.unsetf(std::cout.scientific);
	std::cout << (std::cout.flags() & std::cout.scientific) << "\n";

	std::cout << stringtodd("  -00123.456000e-5", 0) << "\n";
	std::cout << stringtodd("  -00123.456000e-5", -1) << "\n";
	std::cout << stringtodd("  -00123.456000e-5", 1) << "\n";
	std::cout << stringtodd("1.00000000000000000000000000000000000000000000001", -1) << "\n";
	std::cout << stringtodd("1.00000000000000000000000000000000000000000000001", 1) << "\n";
	std::cout << stringtodd("1.00000000000000000000000000000000000000000000001", -1, true) << "\n";
	std::cout << stringtodd("1.00000000000000000000000000000000000000000000001", 1, true) << "\n";
}
