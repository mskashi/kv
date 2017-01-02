/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDOUBLE_HPP
#define RDOUBLE_HPP

#ifdef KV_NOHWROUND
#include <kv/rdouble-nohwround.hpp>
#else
#include <kv/rdouble-hwround.hpp>
#endif

#endif // RDOUBLE_HPP
