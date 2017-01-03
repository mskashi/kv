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

int main()
{
	std::cout << "all test should be 1.\n";
	std::cout << "add_up: " << check_add_up() << "\n";
	std::cout << "add_down: " << check_add_down() << "\n";
	std::cout << "sub_up: " << check_sub_up() << "\n";
	std::cout << "sub_down: " << check_sub_down() << "\n";
	std::cout << "mul_up: " << check_mul_up() << "\n";
	std::cout << "mul_down: " << check_mul_down() << "\n";
	std::cout << "div_up: " << check_div_up() << "\n";
	std::cout << "div_down: " << check_div_down() << "\n";
	std::cout << "sqrt_up: " << check_sqrt_up() << "\n";
	std::cout << "sqrt_down: " << check_sqrt_down() << "\n";
}
