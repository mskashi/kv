/*
 * Copyright (c) 2026 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef FP80_HPP
#define FP80_HPP

/*
 * Determine whether Intel extended 80-bit floating-point is available;
 * if so, define the macro KV_HAVE_FP80
 * and define kv::fp80 as an alias for long double.
 */

#if defined(_MSC_VER)
  /* MSVC is unsupported: long double is the same as double. */
#elif (defined(__i386__) || defined(__x86_64__)) && \
      defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ == 64) && \
      defined(__LDBL_MAX_EXP__)  && (__LDBL_MAX_EXP__  == 16384) && \
      defined(__LDBL_MIN_EXP__)  && (__LDBL_MIN_EXP__  == -16381)
  /* Intel extended 80-bit floating-point is available. */
#define KV_HAVE_FP80 1
#endif

#ifdef KV_HAVE_FP80
namespace kv {
using fp80 = long double;
}
#endif

#endif // FP80_HPP
