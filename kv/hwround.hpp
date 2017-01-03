/*
 * Copyright (c) 2013-2016 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef HWROUND_HPP
#define HWROUND_HPP

/*
 * Hardware Rounding Mode Change
 */

#if defined(KV_FASTROUND)
  #include <emmintrin.h>
#else
  #if defined(_MSC_VER)
    #include <float.h>
  #else
    #include <fenv.h>
  #endif
#endif

namespace kv {

#if defined(KV_FASTROUND)

// Fast implementation for 64bit.
// Assuming that all floating point operation use SSE2. DO NOT USE FPU.
// Assuming that all bits of MXCSR control register except for rounding control
//  is not changed throughout the program execution.

// hwround_regs will be initialized before main()
static struct hwround_reg {
	unsigned long int nearest;
	unsigned long int down;
	unsigned long int up;
	unsigned long int chop;

	hwround_reg() {
		unsigned long int reg = _mm_getcsr();
		nearest = (reg & ~0x00006000) | 0x00000000;
		down = (reg & ~0x00006000) | 0x00002000;
		up =   (reg & ~0x00006000) | 0x00004000;
		chop = (reg & ~0x00006000) | 0x00006000;
	}
} hwround_regs;

struct hwround {
	static void roundnear() {
		_mm_setcsr(hwround_regs.nearest);
	}

	static void rounddown() {
		_mm_setcsr(hwround_regs.down);
	}

	static void roundup() {
		_mm_setcsr(hwround_regs.up);
	}

	static void roundchop() {
		_mm_setcsr(hwround_regs.chop);
	}
};

#else
#if defined(_MSC_VER)

// for Visual C++ (Visual C++ do not have fenv.h)

struct hwround {
	static unsigned int *current_word() {
		static unsigned int cw = 0;
		return &cw;
	}

	static void roundnear() {
		_controlfp_s(current_word(), _RC_NEAR, _MCW_RC);
	}

	static void rounddown() {
		_controlfp_s(current_word(), _RC_DOWN, _MCW_RC);
	}

	static void roundup() {
		_controlfp_s(current_word(), _RC_UP, _MCW_RC);
	}

	static void roundchop() {
		_controlfp_s(current_word(), _RC_UP, _MCW_RC);
	}
};

#else

// for complilers which have fenv.h

struct hwround {
	static void roundnear() {
		fesetround(FE_TONEAREST);
	}

	static void rounddown() {
		fesetround(FE_DOWNWARD);
	}

	static void roundup() {
		fesetround(FE_UPWARD);
	}

	static void roundchop() {
		fesetround(FE_TOWARDZERO);
	}
};

#endif
#endif

} // namespace kv

#endif // HWROUND_HPP
