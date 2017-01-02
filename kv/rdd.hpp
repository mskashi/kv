/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_HPP
#define RDD_HPP

#ifdef KV_NOHWROUND
#include <kv/rdd-nohwround.hpp>
#else
#include <kv/rdd-hwround.hpp>
#endif

#endif // RDD_HPP
