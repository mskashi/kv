#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <cmath>

bool check_add_up()
{
	bool result = true;
	double c1 = pow(2., -54);
	double c2 = 1. + pow(2., -52);
	double c3 = pow(2., -55);
	double c4 = -1. + pow(2., -53);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::add_up(1., c1) == c2);
	result = result && (kv::rop<double>::add_up(-1., c3) == c4);
	kv::rop<double>::end();
	return result;
}

bool check_add_down()
{
	bool result = true;
	double c1 = -pow(2., -55);
	double c2 = 1. - pow(2., -53);
	double c3 = -pow(2., -54);
	double c4 = -1. - pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::add_down(1., c1) == c2);
	result = result && (kv::rop<double>::add_down(-1., c3) == c4);
	kv::rop<double>::end();
	return result;
}

bool check_sub_up()
{
	bool result = true;
	double c1 = -pow(2., -54);
	double c2 = 1. + pow(2., -52);
	double c3 = -pow(2., -55);
	double c4 = -1. + pow(2., -53);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::sub_up(1., c1) == c2);
	result = result && (kv::rop<double>::sub_up(-1., c3) == c4);
	kv::rop<double>::end();
	return result;
}

bool check_sub_down()
{
	bool result = true;
	double c1 = pow(2., -55);
	double c2 = 1. - pow(2., -53);
	double c3 = pow(2., -54);
	double c4 = -1. - pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::sub_down(1., c1) == c2);
	result = result && (kv::rop<double>::sub_down(-1., c3) == c4);
	kv::rop<double>::end();
	return result;
}

bool check_mul_up()
{
	bool result = true;
	double c1 = 1. + pow(2., -52);
	double c2 = 1. + 3. * pow(2., -52);
	double c3 = -1. - 2. * pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::mul_up(c1, c1) == c2);
	result = result && (kv::rop<double>::mul_up(c1, -c1) == c3);
	kv::rop<double>::end();
	return result;
}

bool check_mul_down()
{
	bool result = true;
	double c1 = 1. + pow(2., -52);
	double c2 = 1. + 2. * pow(2., -52);
	double c3 = -1. - 3. * pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::mul_down(c1, c1) == c2);
	result = result && (kv::rop<double>::mul_down(c1, -c1) == c3);
	kv::rop<double>::end();
	return result;
}

bool check_div_up()
{
	bool result = true;
	double c1 = 1. - pow(2., -52);
	double c2 = 1. + 2. * pow(2., -52);
	double c3 = -1. - pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::div_up(1., c1) == c2);
	result = result && (kv::rop<double>::div_up(-1, c1) == c3);
	kv::rop<double>::end();
	return result;
}

bool check_div_down()
{
	bool result = true;
	double c1 = 1. - pow(2., -52);
	double c2 = 1. + pow(2., -52);
	double c3 = -1. - 2. * pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::div_down(1., c1) == c2);
	result = result && (kv::rop<double>::div_down(-1, c1) == c3);
	kv::rop<double>::end();
	return result;
}

bool check_sqrt_up()
{
	bool result = true;
	double c1 = 1. + 3. * pow(2., -52);
	double c2 = 1. + 2. * pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::sqrt_up(c1) == c2);
	kv::rop<double>::end();
	return result;
}

bool check_sqrt_down()
{
	bool result = true;
	double c1 = 1. + 3. * pow(2., -52);
	double c2 = 1. +  pow(2., -52);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::sqrt_down(c1) == c2);
	kv::rop<double>::end();
	return result;
}

bool check_flush_to_zero()
{
	bool result = true;
	double c1 = pow(2., -537);

	kv::rop<double>::begin();
	result = result && (kv::rop<double>::div_up(kv::rop<double>::mul_up(c1, c1), c1) == c1);
	kv::rop<double>::end();
	return result;
}

void check(bool (*f)(), bool& result, const char *str) {
	bool tmp = f();
	result = result && tmp;
	if (tmp) {
		std::cout << str << ": OK\n";
	} else {
		std::cout << str << ": NG\n";
	}
}

int main()
{
	bool result = true;
	check(check_add_up, result, "add_up");
	check(check_add_down, result, "add_down");
	check(check_sub_up, result, "sub_up");
	check(check_sub_down, result, "sub_down");
	check(check_mul_up, result, "mul_up");
	check(check_mul_down, result, "mul_down");
	check(check_div_up, result, "div_up");
	check(check_div_down, result, "div_down");
	check(check_sqrt_up, result, "sqrt_up");
	check(check_sqrt_down, result, "sqrt_down");
	check(check_flush_to_zero, result, "no flush_to_zero");
	if (result) {
		std::cout << "\nall tests passed\n";
	} else {
		std::cout << "\nsome tests failed\n";
	}
}
