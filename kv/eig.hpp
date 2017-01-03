/*
 * eig.hpp:
 * Computing all matrix eigenvalues and all eigenvectors A*V=V*D
 *
 * written	Jun. 4, 2015	A. Takayasu
 * modified	by Masahide Kashiwagi
 * modified Oct. 11, 2015 A. Takayasu
 */

#ifndef EIG_HPP
#define EIG_HPP

#include <iostream>
#include <cmath>
#include <limits>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/complex.hpp>
#include <kv/vleq.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

// Eigenvalue computation using QR method for non-symmetric matrix

namespace kv {

namespace ub = boost::numeric::ublas;

template <class T> bool house(const ub::matrix<T>& x, ub::matrix<T>& v, T& beta)
{
	// function [v,beta] = house(x)
	// n = length(x);
	int n = x.size1();
	int i;
	T sigma, mu;
	v.resize(n,1);
	v(0,0) = 1.;
	sigma = 0.;
	for (i=1; i<n; i++){
		sigma += x(i,0)*x(i,0); // sigma = x(2:n)'*x(2:n);
		v(i,0) = x(i,0); // v = [1;x(2:n)];
	}
	//
	// if sigma == 0
	//   beta = 0;
	// else
	//   mu = sqrt(x(1)^2+sigma);
	//   if x(1)<=0
	//     v(1) = x(1)-mu;
	//   else
	//     v(1) = -sigma/(x(1)+mu);
	//   end
	//   beta = 2*v(1)^2/(sigma+v(1)^2);
	//   v = v/v(1);
	// end
	//
	if (sigma == 0.){
		beta = 0.;
	} else {
		mu = sqrt(x(0,0) * x(0,0) + sigma);
		if (x(0,0) <= 0.){
			v(0,0) = x(0,0) - mu;
		} else {
			v(0,0) = -sigma / (x(0,0) + mu);
		}
		beta = 2 * v(0,0)*v(0,0) / (sigma + v(0,0)*v(0,0));
		v = v / v(0,0);
	}
	return true;
}

template <class T> bool hess(const ub::matrix<T>& A, ub::matrix<T>& Q, ub::matrix<T>& H)
{
	int i,k=0;
	T beta;
	ub::vector<T> vec_tmp;

	// n = size(A,2);
	int n = A.size1();
	if (n != A.size2()) return false;// Square matrix only

	// Q = eye(n);
	ub::matrix<T> x, v, mat_H;
	Q = ub::identity_matrix<T>(n);
	H = A;

	for (k = 0; k < n-2; k++) {
		// [v,beta] = house(H(k+1:n,k));
		x=ub::project(H,ub::range(k+1,n),ub::range(k,k+1));
		house(x,v,beta);
		//   mat_H = eye(n-k) - beta*v*(v');
		mat_H = ub::identity_matrix<T>(n-k-1) - beta*prod(v,trans(v));
		//   H(k+1:n,k) = [mat_H(1,:)*A(k+1:n,k);zeros(n-k-1,1)];
		vec_tmp  = prod(row(mat_H,0),ub::project(H,ub::range(k+1,n),ub::range(k,k+1)));
		H(k+1,k) = vec_tmp(0);
		for (i=k+2;i<n;i++){
			H(i,k)= (T) 0;
		}
		//   H(k+1:n,k+1:n) = mat_H*H(k+1:n,k+1:n);
		ub::project(H,ub::range(k+1,n),ub::range(k+1,n)) = prod(mat_H,ub::project(H,ub::range(k+1,n),ub::range(k+1,n)));
		//   H(1:n,k+1:n) = H(1:n,k+1:n)*mat_H;
		ub::project(H,ub::range(0,n),ub::range(k+1,n)) = prod(ub::project(H,ub::range(0,n),ub::range(k+1,n)),mat_H);
		//   Q(:,k+1:n)=Q(:,k+1:n)*mat_H;
		ub::project(Q,ub::range(0,n),ub::range(k+1,n)) = prod(ub::project(Q,ub::range(0,n),ub::range(k+1,n)),mat_H);
	}
	return true;
}

template <class T> bool francisQR(ub::matrix<T>& Q, ub::matrix<T>& H)
{
	// double tol=1e-15; // Little bit strong!!!
	T tol = std::numeric_limits<T>::epsilon();
	int n = H.size1();
	// % p indicates the ‘active’ matrix size
	int p = n, q, r;
	T s,t,x,y,z,beta;
	ub::matrix<T> v, mat_tmp, mat_tmp2, mat_T;
	mat_tmp.resize(3,1);
	mat_tmp2.resize(2,1);

	while (p > 2) {
		q = p-1;
		s = H(q-1,q-1) + H(p-1,p-1);
		t = H(q-1,q-1)*H(p-1,p-1) - H(q-1,p-1)*H(p-1,q-1);
		// % compute first 3 elements of first column of M
		x = H(0,0)*H(0,0)+H(0,1)*H(1,0)-s*H(0,0)+t;
		y = H(1,0)*(H(0,0)+H(1,1)-s);
		z = H(1,0)*H(2,1);
		for (int k = 0; k < p-2; k++) {
			mat_tmp(0,0) = x;
			mat_tmp(1,0) = y;
			mat_tmp(2,0) = z;
			// [v,beta] = house([x,y,z].');
			house(mat_tmp,v,beta);
			r = fmax(1,k); // r = max(1,k); Need math.h???
			// T = eye(3) - beta*v*(v');
			mat_T = ub::identity_matrix<T>(3) - beta*prod(v,trans(v));
			// H(k+1:k+3,r:n) = T*H(k+1:k+3,r:n);
			ub::project(H,ub::range(k,k+3),ub::range(r-1,n)) = prod(mat_T,ub::project(H,ub::range(k,k+3),ub::range(r-1,n)));
			r = fmin(k+4,p);
			// H(1:r,k+1:k+3) = H(1:r,k+1:k+3)*T;
			ub::project(H,ub::range(0,r),ub::range(k,k+3)) = prod(ub::project(H,ub::range(0,r),ub::range(k,k+3)),mat_T);
			// Q(:,k+1:k+3) = Q(:,k+1:k+3)*T;
			ub::project(Q,ub::range(0,n),ub::range(k,k+3)) = prod(ub::project(Q,ub::range(0,n),ub::range(k,k+3)),mat_T);
			// x = H(k+2,k+1);
			x = H(k+1,k);
			// y = H(k+3,k+1);
			y = H(k+2,k);
			// if k<p-3, z=H(k+4,k+1);
			if (k<p-3) {
				z = H(k+3,k);
			}
		}
		mat_tmp2(0,0) = x;
		mat_tmp2(1,0) = y;
		// [v,beta] = house([x,y]');
		house(mat_tmp2,v,beta);
		// T = eye(2) - beta*v*(v');
		mat_T = ub::identity_matrix<T>(2) - beta*prod(v,trans(v));
		// H(q:p,p-2:n) = T'*H(q:p,p-2:n);
		ub::project(H,ub::range(q-1,p),ub::range(p-3,n)) = prod(mat_T,ub::project(H,ub::range(q-1,p),ub::range(p-3,n)));
		// H(1:p,p-1:p) = H(1:p,p-1:p)*T;
		ub::project(H,ub::range(0,p),ub::range(p-2,p)) = prod(ub::project(H,ub::range(0,p),ub::range(p-2,p)),mat_T);
		// Q(:,q:p) = Q(:,q:p)*T;
		ub::project(Q,ub::range(0,n),ub::range(q-1,p)) = prod(ub::project(Q,ub::range(0,n),ub::range(q-1,p)),mat_T);
		// check for convergence
		// if abs(H(p,q)) < tol*(abs(H(q,q))+abs(H(p,p)))
		// if (fabs(H(p-1,q-1)) < tol*(fabs(H(q-1,q-1) + fabs(H(p-1,p-1)))))
		using std::abs;
		if (abs(H(p-1,q-1)) < tol*(abs(H(q-1,q-1) + abs(H(p-1,p-1))))) {
			H(p-1,q-1) = (T) 0;
			p--;
		} else if (abs(H(p-2,q-2)) < tol*(abs(H(q-2,q-2)) + fabs(H(q-1,q-1)))){
			H(p-2,q-2) = (T) 0;
			p-=2;
		}
	}
	return true;
}

template <class T> bool lu_factorize_comp(ub::matrix<kv::complex<T> >& A, ub::vector<int>& p)// Numerical recipe in C (ludcmp)
{
  int i,j,k,imax;
  int n = A.size1();

  ub::vector<T> vv;
  vv.resize(n);

  T big, temp;
  kv::complex<T> sum, dum;

  using std::abs;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++) {
      if ((temp = abs(A(i,j))) > big) {
        big = temp;
      }
    }
    if (big == 0.0) {
      std::cout << "Singular matrix in lu_factorize_comp" << std::endl;
      return false;
    }
    vv(i) = 1.0/big;
  }

  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = A(i,j);
      for (k = 0; k < i; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = A(i,j);
      for (k = 0; k < j; k++) {
        sum -= A(i,k)*A(k,j);
      }
      A(i,j) = sum;
      if ((temp=vv(i)*abs(sum)) >= big) {
        big = temp;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = A(imax,k);
        A(imax,k) = A(j,k);
        A(j,k) = dum;
      }
      vv(imax) = vv(j);
    }
    p(j) = imax;
    if (abs(A(j,j)) == 0) {
      std::cout << "Singular matrix in lu_factorize_comp" << std::endl;
      return false;
    }
    if (j != n) {
      dum = 1./(A(j,j));
      for (i = j+1; i < n; i++) {
        A(i,j) *= dum;
      }
    }
  }
  return true;
}

template <class T> bool lu_substitute_comp(const ub::matrix<kv::complex<T> >& A, const ub::vector<int>& p, ub::vector<kv::complex<T> >& b)// Numerical recipe in C (lubksb)
{
  int i, ii=0, ip, j;
  int n = A.size1();
  kv::complex<T> sum;

  using std::abs;

  for (i = 0; i < n; i++) {
    ip = p(i);
    sum = b(ip);
    b(ip) = b(i);
    if (ii==0) {
      for (j = ii; j <= i-1; j++) {
        sum -= A(i,j)*b(j);
      }
    } else if (abs(sum) != 0) {
      ii = i;
    }
    b(i) = sum;
  }
  for (i = n-1; i >= 0; i--) {
      sum = b(i);
      for (j = i+1; j < n; j++) {
        sum -= A(i,j)*b(j);
      }
      b(i) = sum/A(i,i);
  }
  return true;
}

template <class T> bool eig2by2(const ub::matrix<T>& P, const ub::matrix<T>& A, ub::matrix< kv::complex<T> >& V, ub::matrix< kv::complex<T> >& D)
{
	int i, j, k, l;
	T tra, det, x_norm;
	kv::complex<T> tmp, tmp1, am, ap;
	kv::complex<T> b1, b2;
	kv::complex<T> a11, a12, a21, a22, ck;
	int n = A.size1();// n = size(A,2);
	// X = eye(n,n); Orthogonal matrix
	ub::matrix< kv::complex<T> > X;
	ub::vector< kv::complex<T> > x;
	X.resize(n,n);
	for (i = 0; i < n; i++) {
		X(i,i) = 1.;
	}
	D.resize(n,n);// Eigen value matrix (diagonal)

	i=1;
	// Compute eigenvalue
	while (i<=n) {
		if (i!=n) {
			j = i+1;
		}
		if (A(j-1,i-1) == 0 || i==n) {
			// Eigen value is diagonal element
			D(i-1,i-1) = A(i-1,i-1);
			i++;
		} else {
			// If A contains 2 by 2 block on the diagonal,
			// compute a real pair or a complex conjugate pair.
			tra = A(i-1,i-1) + A(j-1,j-1);
			det = A(i-1,i-1)*A(j-1,j-1)-A(i-1,j-1)*A(j-1,i-1);
			tmp = tra*tra-4*det;// Complex!
			tmp1 = sqrt(tmp);
			am = 0.5*(tra + tmp1);
			ap = 0.5*(tra - tmp1);
			using std::abs;
			if (abs(A(i-1,i-1)-am)/abs(am) < 1) {
				D(i-1,i-1) = am;
				D(j-1,j-1) = ap;
			} else {
				D(i-1,i-1) = ap;
				D(j-1,j-1) = am;
			}
			i += 2;
		}
	}
	// std::cout << "D" << D << std::endl;
	// D(0,0).imag() = -D(0,0).imag();
	// std::cout << "conj(D)" << D << std::endl;


	// Compute eigenvector by backward substitution
	bool flag = true;
	for (i = 1; i <= n; i++) {
		j = i;
		if (flag) {
			// Compute jth element of ith eigenvector X(j,i)
			while (j>0) {
				k = fmin(n,j+1);
				l = fmax(1,j-1);
				if (i==j) {
					if (A(k-1,j-1) != 0 && k!=j) {
						// Block diagonal [A(j,j), A(j,j+1); A(j+1,j), A(j+1,j+1)] appears
						using std::abs;
						if (abs(A(j,j)-D(i-1,i-1)) > abs(A(j-1,j))) {
							X(j,i-1) = -A(j,j-1) / (A(j,j)-D(i-1,i-1));
						} else {
							X(j,i-1) = -(A(j-1,j-1)-D(i-1,i-1))/A(j-1,j);
						}
						if (D(j-1,j-1).imag() != -D(j-1,j-1).imag()) {
							// Conjugate pair X(j,i) = 1+0*1i;
							flag = false; //Conjugate pair flag (flag=false means the next eigenvector is conjugate of ith vector)
						}
						j--;
					} else if (A(j-1,l-1) != 0 && l!=j) {
						// Block diagonal [A(j-1,j-1), A(j-1,j); A(j,j-1), A(j,j)] appears X(j,i)=1
						using std::abs;
						if (abs(A(j-2,j-2)-D(i-1,i-1)) > abs(A(j-1,j-2))) {
							X(j-2,i-1) = -A(j-2,j-1)/(A(j-2,j-2)-D(i-1,i-1));
						} else {
							X(j-2,i-1) = -(A(j-1,j-1)-D(i-1,i-1))/A(j-1,j-2);
						}
						if (D(j-1,j-1).imag() != -D(j-1,j-1).imag()) {
							// Conjugate pair X(j,j) = 1+0*1i;
							flag = false;// Conjugate pair flag (flag=false means the next eigenvector is conjugate of ith vector)
						}
						j -= 2;
					} else {
						// A(i,i) is eigen value
						j--;
					}
				} else {
					if (A(j-1,l-1)!=0 && l!=j) {
						// Block diagonal [A(j-1,j-1), A(j-1,j); A(j,j-1), A(j,j)] appears
						// Real pair
						b1=0; b2=0;
						for (k = j+1; k <= fmin(i+1,n); k++) {
							b1 += A(j-2,k-1)*X(k-1,i-1);
							b2 += A(j-1,k-1)*X(k-1,i-1);
						}
						a11 = A(j-2,j-2)-D(i-1,i-1);
						a12 = A(j-2,j-1);
						a21 = A(j-1,j-2);
						a22 = A(j-1,j-1)-D(i-1,i-1);
						ub::matrix<kv::complex<T> > LU(2,2);
						LU(0,0) = a11; LU(0,1) = a12; LU(1,0) = a21; LU(1,1) = a22;
						ub::vector<kv::complex<T> > b(2);
						b(0) = b1; b(1) = b2;
						ub::vector<int> pm(n);
						lu_factorize_comp(LU,pm);
						lu_substitute_comp(LU,pm,b);
						// ck  = -1/(a11*a22-a12*a21);
	          // X(j-2,i-1) = ck*(a22*b1-a12*b2);
	          // X(j-1,i-1) = ck*(-a21*b1+a11*b2);
						X(j-2,i-1) = b(0);
						X(j-1,i-1) = b(1);
						j -= 2;
					} else {
						// A(i,i) is eigen value
						for (k = j+1; k <= fmin(i+1,n); k++) {
							X(j-1,i-1) = X(j-1,i-1) + A(j-1,k-1)*X(k-1,i-1);
						}
						// X(j-1,i-1) = -X(j-1,i-1)/(D(j-1,j-1)-D(i-1,i-1));
						if ((D(j-1,j-1)-D(i-1,i-1)).real()==0 && (D(j-1,j-1)-D(i-1,i-1)).imag()==0 && X(j-1,i-1).real()==0 && X(j-1,i-1).imag()==0) {
							X(j-1,i-1) = 0; // d(j) = d(i)
						} else {
							X(j-1,i-1) = -X(j-1,i-1)/(D(j-1,j-1)-D(i-1,i-1));
						}
						j--;
					}
				}
			}
			x = column(X,i-1);
			x_norm = 0;
			for (k = 0; k < n; k++) {
				using std::abs;
				x_norm += abs(x(k))*abs(x(k));
			}
			x_norm = sqrt(x_norm);
			for (k = 0; k < n; k++) {
				X(k,i-1) /= x_norm;
			}
		} else {
			for (k = 0; k < n; k++) {
				X(k,i-1).real() =  X(k,i-2).real();
				X(k,i-1).imag() =  -X(k,i-2).imag();
			}
			flag = true;
		}
	}
	// V = P*X;
	V = prod(P,X);

	return true;
}

template <class T> bool eig(const ub::matrix<T>& A, ub::matrix< kv::complex<T> >& V, ub::matrix< kv::complex<T> >& D)
{
	int n = A.size1();
	if (n != A.size2()) return false;// Square matrix only

	ub::matrix<T> Q, H;
	Q.resize(n,n);
	H.resize(n,n);

	// std::cout << Q << "\n";
	hess(A,Q,H);
	francisQR(Q,H);
	eig2by2(Q,H,V,D);
	// std::cout << "Residual " << prod(A,Q)-prod(Q,H) << "\n";
	// std::cout << "Residual " << prod(A,V)-prod(V,D) << "\n";
	return true;
}

template <class T> bool veig(const ub::matrix<T>& A, ub::vector< kv::complex< kv::interval<T> > >& v)
{
	int i, j, n=A.size1();
	ub::matrix< kv::complex<T> > V, D;
	ub::matrix< kv::complex< kv::interval<T> > > X, C;
	ub::matrix< kv::interval<T> > BA, BC, G;
	// ub::vector< kv::complex< kv::interval<T> > > b, x;
	ub::vector< kv::interval<T> > d, err;

	eig(A,V,D);
	// std::cout << "Residual " << prod(a,V)-prod(V,D) << "\n";
	// C = A*intval(X);
	X.resize(n,n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			X(i,j) = V(i,j);
		}
	}
	C = prod(A,X);

	// G = verifylss(X,C);
	G.resize(2*n,n);
	BA.resize(2*n,2*n);
	BC.resize(2*n,n);

	// X  = A + Bi;
	// BA = [A, -B; B, A]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			BA(i,j) = X(i,j).real();
			BA(i+n,j) = X(i,j).imag();
			BA(i,j+n) = -X(i,j).imag();
			BA(i+n,j+n) = X(i,j).real();
		}
	}

	// C = P + Qi
	// BC = [P;Q]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			BC(i,j) = C(i,j).real();
			BC(i+n,j) = C(i,j).imag();
		}
	}
	kv::vleq(BA,BC,G);

	// mid_G = mid(G);

	d.resize(2*n);
	err.resize(2*n);

	for (i = 0; i < n; i++) {
		d(i) = mid(G(i,i));
		d(i+n) = mid(G(i+n,i));
		G(i,i) -= d(i);
		G(i+n,i) -= d(i+n);
	}

	// G = G-mid_G;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			err(i) += G(i,j);
			err(i+n) += G(i+n,j);
		}
	}

	v.resize(n);
	d += err;
	for (i = 0; i < n; i++) {
		v(i) = kv::complex< kv::interval<T> >(d(i),d(i+n));
	}
	return true;
}

template <class T> bool veig(const ub::matrix< kv::interval<T> >& A, ub::vector< kv::complex< kv::interval<T> > >& v)
{
	int i, j, n=A.size1();
	ub::matrix< kv::complex<T> > V, D;
	ub::matrix< kv::complex< kv::interval<T> > > X, C;
	ub::matrix< kv::interval<T> > BA, BC, G;
	// ub::vector< kv::complex< kv::interval<T> > > b, x;
	ub::vector< kv::interval<T> > d, err;

	eig(mid(A),V,D);
	// std::cout << "Residual " << prod(a,V)-prod(V,D) << "\n";
	// C = A*intval(X);
	X.resize(n,n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			X(i,j) = V(i,j);
		}
	}
	C = prod(A,X);

	// G = verifylss(X,C);
	G.resize(2*n,n);
	BA.resize(2*n,2*n);
	BC.resize(2*n,n);

	// X  = A + Bi;
	// BA = [A, -B; B, A]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			BA(i,j) = X(i,j).real();
			BA(i+n,j) = X(i,j).imag();
			BA(i,j+n) = -X(i,j).imag();
			BA(i+n,j+n) = X(i,j).real();
		}
	}

	// C = P + Qi
	// BC = [P;Q]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			BC(i,j) = C(i,j).real();
			BC(i+n,j) = C(i,j).imag();
		}
	}
	kv::vleq(BA,BC,G);

	// mid_G = mid(G);

	d.resize(2*n);
	err.resize(2*n);

	for (i = 0; i < n; i++) {
		d(i) = mid(G(i,i));
		d(i+n) = mid(G(i+n,i));
		G(i,i) -= d(i);
		G(i+n,i) -= d(i+n);
	}

	// G = G-mid_G;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			err(i) += G(i,j);
			err(i+n) += G(i+n,j);
		}
	}

	v.resize(n);
	d += err;
	for (i = 0; i < n; i++) {
		v(i) = kv::complex< kv::interval<T> >(d(i),d(i+n));
	}
	return true;
}

template <class T> bool invert_comp(const ub::matrix<kv::complex<T> >& A, ub::matrix<kv::complex<T> >& R)
{
  int n = A.size1();
  ub::matrix< kv::complex<T> > LU=A;

  ub::vector< kv::complex<T> > b(n);
  ub::vector<int> pm(n);

  R.resize(n,n);

  lu_factorize_comp(LU,pm);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j==i) b(j) = 1.0;
      else b(j) = 0.0;
    }
    // std::cout << "b" << b << std::endl;
    lu_substitute_comp(LU,pm,b);
    ub::column(R,i) = b;
  }
  return true;
}

} // namespace kv

#endif // EIG_HPP
