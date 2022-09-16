/*
 * Copyright (c) 2013-2022 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef HWROUND_HPP
#define HWROUND_HPP

/*
 * Hardware Rounding Mode Change
 */

#if defined(KV_FASTROUND)

#include <cstdint>

#if defined(_WIN32) || defined(_WIN64) || defined(__i386__) || defined(__x86_64__)
#include <emmintrin.h>
#endif

#else // KV_FASTROUND

#if defined(_MSC_VER)
#include <float.h>
#else
#include <fenv.h>
#endif

#endif // KV_FASTROUND


namespace kv {

#if defined(KV_FASTROUND)
  #if defined(_WIN32) || defined(_WIN64) || defined(__i386__) || defined(__x86_64__)

// Fast implementation of rounding-mode-change for SSE2
// Assuming that all floating point operation use SSE2. DO NOT USE FPU.
// Assuming that all bits of MXCSR control register except for rounding control
//  is not changed throughout the program execution.

// hwround_regs will be initialized before main()
static struct hwround_reg {
	uint32_t nearest;
	uint32_t down;
	uint32_t up;
	uint32_t chop;

	hwround_reg() {
		uint32_t reg = _mm_getcsr();
		nearest = (reg & ~(3UL << 13)) | (0UL << 13);
		down = (reg & ~(3UL << 13)) | (1UL << 13);
		up =   (reg & ~(3UL << 13)) | (2UL << 13);
		chop = (reg & ~(3UL << 13)) | (3UL << 13);
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

#endif

#if defined(__aarch64__)

// Assuming that all bits of FPCR register except for rounding control
//  is not changed throughout the program execution.

// hwround_regs will be initialized before main()
static struct hwround_reg {
	uint64_t nearest;
	uint64_t down;
	uint64_t up;
	uint64_t chop;

	hwround_reg() {
		uint64_t reg;
		asm volatile ("mrs %0, fpcr" : "=r" (reg));
		nearest = (reg & ~(3ULL << 22)) | (0ULL << 22);
		down = (reg & ~(3ULL << 22)) | (2ULL << 22);
		up = (reg & ~(3ULL << 22)) | (1ULL << 22);
		chop = (reg & ~(3ULL << 22)) | (3ULL << 22);
	}
} hwround_regs;

struct hwround {
	static void roundnear() {
		asm volatile ("msr fpcr, %0" : : "r" (hwround_regs.nearest));
	}

	static void rounddown() {
		asm volatile ("msr fpcr, %0" : : "r" (hwround_regs.down));
	}

	static void roundup() {
		asm volatile ("msr fpcr, %0" : : "r" (hwround_regs.up));
	}

	static void roundchop() {
		asm volatile ("msr fpcr, %0" : : "r" (hwround_regs.chop));
	}
};
#endif // __aarch64__

#if defined(__arm__)

// Assuming that all bits of FPSCR register except for rounding control
//  is not changed throughout the program execution.

// hwround_regs will be initialized before main()
static struct hwround_reg {
	uint32_t nearest;
	uint32_t down;
	uint32_t up;
	uint32_t chop;

	hwround_reg() {
		uint32_t reg;
		asm volatile ("vmrs %0, fpscr" : "=r" (reg));
		nearest = (reg & ~(3UL << 22)) | (0UL << 22);
		down = (reg & ~(3UL << 22)) | (2UL << 22);
		up = (reg & ~(3UL << 22)) | (1UL << 22);
		chop = (reg & ~(3UL << 22)) | (3UL << 22);
	}
} hwround_regs;

struct hwround {
	static void roundnear() {
		asm volatile ("vmsr fpscr, %0" : : "r" (hwround_regs.nearest));
	}

	static void rounddown() {
		asm volatile ("vmsr fpscr, %0" : : "r" (hwround_regs.down));
	}

	static void roundup() {
		asm volatile ("vmsr fpscr, %0" : : "r" (hwround_regs.up));
	}

	static void roundchop() {
		asm volatile ("vmsr fpscr, %0" : : "r" (hwround_regs.chop));
	}
};
#endif // __arm__

#else // KV_FASTROUND

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

#else // _MSC_VER

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

#endif // _MSC_VER

#endif // KV_FASTROUND

} // namespace kv

#endif // HWROUND_HPP
