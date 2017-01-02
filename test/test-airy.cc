#include <kv/airy.hpp>

int main()
{
	std::cout.precision(17);
	int i;

	for (i=-10; i<=10; i++) {
		std::cout << "ai: " << i << ": " << kv::airy_Ai(kv::interval<double>(i)) << "\n";;
	}

	for (i=-10; i<=10; i++) {
		std::cout << "bi: " << i << ": " << kv::airy_Bi(kv::interval<double>(i)) << "\n";;
	}

	for (i=-10; i<=10; i++) {
		std::cout << "aid: " << i << ": " << kv::airy_Ai_d(kv::interval<double>(i)) << "\n";;
	}

	for (i=-10; i<=10; i++) {
		std::cout << "bid: " << i << ": " << kv::airy_Bi_d(kv::interval<double>(i)) << "\n";;
	}
}
