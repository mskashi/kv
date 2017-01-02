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

struct hwround {
	static const unsigned long int& default_register() {
		const static unsigned long int reg = _mm_getcsr();
		return reg;
	}

	static void roundnear() {
		const static unsigned long int reg = (default_register() & ~0x00006000) | 0x00000000;
		_mm_setcsr(reg);
	}

	static void rounddown() {
		const static unsigned long int reg = (default_register() & ~0x00006000) | 0x00002000;
		_mm_setcsr(reg);
	}

	static void roundup() {
		const static unsigned long int reg = (default_register() & ~0x00006000) | 0x00004000;
		_mm_setcsr(reg);
	}

	static void roundchop() {
		const static unsigned long int reg = (default_register() & ~0x00006000) | 0x00006000;
		_mm_setcsr(reg);
	}
};

#else
#if defined(_MSC_VER)

// for Visual C++ (Visual C++ do not have fenv.h)

struct hwround {
	static void roundnear() {
		_controlfp(_RC_NEAR, _MCW_RC);
	}

	static void rounddown() {
		_controlfp(_RC_DOWN, _MCW_RC);
	}

	static void roundup() {
		_controlfp(_RC_UP, _MCW_RC);
	}

	static void roundchop() {
		_controlfp(_RC_UP, _MCW_RC);
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
