#include "tntsupp.h"

using namespace std;
using namespace TNT;

template <class T>
T norm2(const Vector<T> &A) {
  return dot_prod(A, A);
}

template <class T>
T norm2(const Fortran_Matrix<T> &A) {
  Vector<T> v(A.num_cols() * A.num_rows(), A.begin());
  return dot_prod(v, v);
}

// trasform rect cord to hexa cord;
template <class T>
Vector<T> hexa2rect(const Vector<T> &x) {
  Vector<T> ans(x.size());
  if (((int) x(2) % 2) == 0) ans(1) = x(1);
  else ans(1) = x(1) + 0.5;
  ans(2) = sqrt(3.0) / 2 * x(2);
  return ans;
}

template <class T> 
Vector<T> hexa2rect(const Fortran_Matrix<T> &x) {
  // x is 1 row
  Vector<T> ans(x.num_cols());
  if (((int) x(1,2) % 2) == 0) ans(1) = x(1,1);
  else ans(1) = x(1,1) + 0.5;
  ans(2) = sqrt(3.0) / 2 * x(1,2);
  return ans;
}

template <class T> 
Vector<T> hexa2rect(const Region2D<Fortran_Matrix<T> > &x) {
  // x is 1 row
  Vector<T> ans(x.num_cols());
  if (((int) x(1,2) % 2) == 0) ans(1) = x(1,1);
  else ans(1) = x(1,1) + 0.5;
  ans(2) = sqrt(3.0) / 2 * x(1,2);
  return ans;
}

template <class T> 
Vector<T> hexa2rect(const_Region2D<Fortran_Matrix<T> > &x) {
  // x is 1 row
  Vector<T> ans(x.num_cols());
  if (((int) x(1,2) % 2) == 0) ans(1) = x(1,1);
  else ans(1) = x(1,1) + 0.5;
  ans(2) = sqrt(3.0) / 2 * x(1,2);
  return ans;
}

int cord2index(double x, double y, int xdim) {
  return (int) (x + y * xdim) + 1;
}
