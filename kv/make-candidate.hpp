/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef MAKE_CANDIDATE_HPP
#define MAKE_CANDIDATE_HPP

#include <boost/numeric/ublas/vector.hpp>

/*
 * calculate radius of candidate set from newton step vector
 */

namespace kv {

template <class T> void make_candidate(boost::numeric::ublas::vector<T>& x)
{
	int n = x.size();
	int i;


#ifdef NEWTON_TWICE
	T m(0.);
	for (i=0; i<n; i++) {
		if (x(i) > m) m = x(i);
	}
	for (i=0; i<n; i++) x(i) = 2. * m;
#else
	T av(0.);

	for (i=0; i<n; i++) av += x(i);

	av /= n;

	for (i=0; i<n; i++) x(i) += av;
#endif
}

} // namespace kv

#endif // MAKE_CANDIDATE_HPP
