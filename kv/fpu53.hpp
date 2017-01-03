/*
 * Copyright (c) 2013-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef FPU53_HPP
#define FPU53_HPP

#if defined(_WIN32) || defined(_WIN64)
#include <float.h>
#endif


namespace kv {

// set FPU's precision to 53bit.
// You *must* run fpu53() if you use dd.hpp on FPU of Intel CPU.
inline void fpu53() {
	#if defined(_WIN32) || defined(_WIN64)
	unsigned int cw = 0;
	_controlfp_s(&cw, _PC_53, _MCW_PC);
	#endif
	#if defined(__i386__) || defined(__x86_64__)
	unsigned short int reg;
	__asm __volatile__ ("fnstcw %0" : "=m"(reg));
	reg &= ~0x0300;
	reg |= 0x0200;
	__asm __volatile__ ("fldcw %0" : : "m"(reg));
	#endif
}

// hack for calling fpu53 before main()
namespace {
	struct dummy_caller {
		dummy_caller() {
			fpu53();
		}
	} dummy;
}

} // namespace kv

#endif // FPU53_HPP
