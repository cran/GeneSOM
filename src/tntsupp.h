#include <iostream>
#include "tnt/subscript.h"
#include "tnt/tnt.h"
#include "tnt/vec.h"
#include "tnt/fmat.h"
#include "tnt/region1d.h"
#include "tnt/region2d.h"
#include "tnt/transv.h"
#include "tnt/lu.h"

using namespace TNT;

typedef Vector<double> DVector;
typedef Vector<int> IVector;
typedef Fortran_Matrix<double> DMatrix;


//typedef Region1D<DVector> subDVector;
//typedef Region2D<DMatrix> subDMatrix;

//matrix operation on Region2D;

template <class T>
Fortran_Matrix<T> operator+(const_Region2D<Fortran_Matrix<T> > &A, 
			    const_Region2D<Fortran_Matrix<T> > &B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M==B.num_rows());
  assert(N==B.num_cols());

  Fortran_Matrix<T> ans(M, N);
  for (Subscript i = 1; i <= M; i++) 
    for (Subscript j = 1; j <= N; j++)
      ans(i, j) = A(i, j) + B(i, j);
  return ans;
}

template <class T>
Fortran_Matrix<T> operator+(const Region2D<Fortran_Matrix<T> > &A, 
			    const Region2D<Fortran_Matrix<T> > &B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M==B.num_rows());
  assert(N==B.num_cols());

  Fortran_Matrix<T> ans(M, N);
  for (Subscript i = 1; i <= M; i++) 
    for (Subscript j = 1; j <= N; j++)
      ans(i, j) = A(i, j) + B(i, j);
  return ans;
}

template <class T>
Fortran_Matrix<T> operator-(const_Region2D<Fortran_Matrix<T> > &A, 
			    const_Region2D<Fortran_Matrix<T> > &B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M==B.num_rows());
  assert(N==B.num_cols());

  Fortran_Matrix<T> ans(M, N);
  for (Subscript i = 1; i <= M; i++) 
    for (Subscript j = 1; j <= N; j++)
      ans(i, j) = A(i, j) - B(i, j);
  return ans;
}
template <class T>
Fortran_Matrix<T> operator-(const Region2D<Fortran_Matrix<T> > &A, 
			    const Region2D<Fortran_Matrix<T> > &B) {
  Subscript M = A.num_rows();
  Subscript N = A.num_cols();

  assert(M==B.num_rows());
  assert(N==B.num_cols());

  Fortran_Matrix<T> ans(M, N);
  for (Subscript i = 1; i <= M; i++) 
    for (Subscript j = 1; j <= N; j++)
      ans(i, j) = A(i, j) - B(i, j);
  return ans;
}

template <class T>
Fortran_Matrix<T> mult_element(const Region2D<Fortran_Matrix<T> > &A, 
			       const Region2D<Fortran_Matrix<T> > &B)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Fortran_Matrix<T> tmp(M,N);
    Subscript i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) * B(i,j);

    return tmp;
}


template <class T>
Fortran_Matrix<T> transpose(const Region2D<Fortran_Matrix<T> > &A)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    Fortran_Matrix<T> S(N,M);
    Subscript i, j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            S(j,i) = A(i,j);

    return S;
}

    
template <class T>
inline Fortran_Matrix<T> matmult(const Region2D<Fortran_Matrix<T> > &A, 
				 const Region2D<Fortran_Matrix<T> > &B)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == B.num_rows());
#endif

    Subscript M = A.num_rows();
    Subscript N = A.num_cols();
    Subscript K = B.num_cols();

    Fortran_Matrix<T> tmp(M,K);
    T sum;

    for (Subscript i=1; i<=M; i++)
    for (Subscript k=1; k<=K; k++)
    {
        sum = 0;
        for (Subscript j=1; j<=N; j++)
            sum = sum +  A(i,j) * B(j,k);

        tmp(i,k) = sum; 
    }

    return tmp;
}

template <class T>
inline Fortran_Matrix<T> operator*(const Region2D<Fortran_Matrix<T> > &A, 
				   const Region2D<Fortran_Matrix<T> > &B)
{
    return matmult(A,B);
}

template <class T>
inline int matmult(Fortran_Matrix<T>& C, const Region2D<Fortran_Matrix<T> >&A, 
		   const Region2D<Fortran_Matrix<T> >&B)
{

    assert(A.num_cols() == B.num_rows());

    Subscript M = A.num_rows();
    Subscript N = A.num_cols();
    Subscript K = B.num_cols();

    C.newsize(M,K);         // adjust shape of C, if necessary


    T sum; 

    const T* row_i;
    const T* col_k;

    for (Subscript i=1; i<=M; i++)
    {
        for (Subscript k=1; k<=K; k++)
        {
            row_i = &A(i,1);
            col_k = &B(1,k);
            sum = 0;
            for (Subscript j=1; j<=N; j++)
            {
                sum +=  *row_i * *col_k;
                row_i += M;
                col_k ++;
            }
        
            C(i,k) = sum; 
        }

    }

    return 0;
}


template <class T>
Vector<T> matmult(const Region2D<Fortran_Matrix<T> > &A, const Vector<T> &x)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == x.dim());
#endif

    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    Vector<T> tmp(M);
    T sum;

    for (Subscript i=1; i<=M; i++)
    {
        sum = 0;
        for (Subscript j=1; j<=N; j++)
            sum = sum +  A(i,j) * x(j);

        tmp(i) = sum; 
    }

    return tmp;
}

template <class T>
inline Vector<T> operator*(const Region2D<Fortran_Matrix<T> > &A, const Vector<T> &x)
{
    return matmult(A,x);
}

template <class T>
inline Fortran_Matrix<T> operator*(const Region2D<Fortran_Matrix<T> > &A, 
				   const T &x)
{
    Subscript M = A.num_rows();
    Subscript N = A.num_cols();

    Subscript MN = M*N; 

    Fortran_Matrix<T> res(M,N);
    const T* a = A.begin();
    T* t = res.begin();
    T* tend = res.end();

    for (t=res.begin(); t < tend; t++, a++)
        *t = *a * x;

    return res;
} 


//convert Region2D to matrix or vector
template <class T>
Fortran_Matrix<T> asMat(const Region2D<Fortran_Matrix<T> > &A) {
  Subscript m = A.num_rows(), n = A.num_cols();
  Fortran_Matrix<T> ans(m, n);
  for (Subscript i = 1; i <= m; i++)
    for (Subscript j = 1; j <= n; j++)
      ans(i, j) = A(i, j);
  return ans;
}

template <class T>
Fortran_Matrix<T> asMat(const_Region2D<Fortran_Matrix<T> > &A) {
  Subscript m = A.num_rows(), n = A.num_cols();
  Fortran_Matrix<T> ans(m, n);
  for (Subscript i = 1; i <= m; i++)
    for (Subscript j = 1; j <= n; j++)
      ans(i, j) = A(i, j);
  return ans;
}

template <class T>
Vector<T> asVec(const Region2D<Fortran_Matrix<T> > &A) {
  // A is 1 row or 1 col
  Subscript m = A.num_rows(), n = A.num_cols();
  if (m == 1) {
    Vector<T> ans(n);
    for (Subscript i = 1; i <= n; i++) ans(i) = A(1,i);
    return ans;
  }
  else {
    Vector<T> ans(m);
    for (Subscript i = 1; i <= m; i++) ans(i) = A(i,1);
    return ans;
  }
}

template <class T>
Vector<T> asVec(const_Region2D<Fortran_Matrix<T> > A) {
  // A is 1 row or 1 col
  Subscript m = A.num_rows(), n = A.num_cols();
  if (m == 1) {
    Vector<T> ans(n);
    for (Subscript i = 1; i <= n; i++) ans(i) = A(1,i);
    return ans;
  }
  else {
    Vector<T> ans(m);
    for (Subscript i = 1; i <= m; i++) ans(i) = A(i,1);
    return ans;
  }
}

//convert vector to matrix
template <class T>
Fortran_Matrix<T> asRowMat(const Vector<T> &v) {
  Subscript n = v.size();
  Fortran_Matrix<T> ans(1,n);
  for (Subscript i = 1; i <= n; i++) ans(1,i) = v(i);
  return ans;
}

template <class T>
Fortran_Matrix<T> asColMat(const Vector<T> &v) {
  Subscript n = v.size();
  Fortran_Matrix<T> ans(n,1);
  for (Subscript i = 1; i <= n; i++) ans(i,1) = v(i);
  return ans;
}


//scalar multiplication 
template <class T>
inline Vector<T> operator*(const Vector<T> &v, const T &x) {
  Subscript m = v.size();
  Vector<T> ans(m);
  for (Subscript i = 1; i <= m; i++)
    ans(i) = v(i) * x;
  return ans;
}

template <class T>
inline Vector<T> operator*(const T &x, const Vector<T> &v) {
  return v * x;
}

template <class T>
inline Fortran_Matrix<T> operator*(const T &x, const Fortran_Matrix<T>  &A) {
  return A * x;
}

// utilities:
Fortran_Matrix<double> ident (int n) {
  Fortran_Matrix<double> ans(n,n);
  for (int i = 1; i <= n; i++) ans(i,i) = 1.0;
  return ans;
}

Fortran_Matrix<double> MatRowCol(const Fortran_Matrix<double> &mat, const Vector<double> &r, const Vector<double> &c) {
  int m = r.size(), n = c.size();
  Fortran_Matrix<double> ans(m,n);
  for (int i = 1; i <= m; i++)
    for (int j = 1; j <= n; j++) 
      ans(i,j) = mat((int) r(i), (int) c(j));
  return ans;
}

template <class T>
Region2D<Fortran_Matrix<T> > MatRow(Fortran_Matrix<T> &x, int i) {
  int n = x.num_cols();
  Index1D J(1,n), I(i,i);
  return x(I,J);
}

template <class T>
Region2D<Fortran_Matrix<T> > MatCol(Fortran_Matrix<T> &x, int i) {
  int m = x.num_rows();
  Index1D I(1,m), J(i,i);
  return x(I,J);
}

template <class T>
Region2D<Fortran_Matrix<T> > MatRows(Fortran_Matrix<T> &x, const Index1D &I) {
  int n = x.num_cols(); Index1D J(1,n);
  return x(I,J);
}

template <class T>
Region2D<Fortran_Matrix<T> > MatCols(Fortran_Matrix<T> &x, const Index1D &J) {
  int m = x.num_rows(); Index1D I(1,m);
  return x(I,J);
}

template <class T>
Region1D<Vector<T> > VecSubs(Vector<T> &x, const Index1D &I) {
  return Region1D<Vector<T> >(x, I);
}

template <class T>
Vector<T> asVec(const Region1D<Vector<T> > &x) {
  Vector<T> ans(x.dim());
  for (int i = 1; i <= ans.size(); i++) ans(i) = x(i);
  return ans;
}


Fortran_Matrix<double> rho2mat(const Vector<double> &rho) {
  int s = rho.size(); // s = n(n-1)/2
  int n = (int) (0.5 * ( 1 + sqrt(1 + 8 * s)));
  Fortran_Matrix<double> fullmat = ident(n); 
  int k = 1;
  for (int i = 1; i <= n - 1; i++)
    for (int j = i + 1; j <= n; j++) {
	fullmat(i, j) = rho(k++);
	fullmat(j, i) = fullmat(i, j);
    }
  return fullmat;
}

//solve(a, b = ident(n))
DMatrix solve(const DMatrix &a, const DMatrix &b) {
  Subscript m = a.dim(1); assert(m == a.dim(2));
  Subscript n = b.dim(1); assert(m == n);
  Subscript l = b.dim(2);
  Vector<Subscript> index(m);
  DMatrix T(a), B(b);
  DMatrix ans(n,l);
  if (LU_factor(T, index) != 0) {
    cout << "LU_factor() failed." << endl; 
    return ans;
  }
  DVector v(m);
  for (int i  = 1; i <= l; i++) {
    v = asVec(MatCol(B,i));
    LU_solve(T, index, v);
    MatCol(ans, i) = asColMat(v);
  }
  return ans;  
}

DVector solve(const DMatrix &A, const DVector &b) {
  DMatrix T(A); Vector<Subscript> index(b.size());
  DVector ans(b);
  if (LU_factor(T, index) !=0) {
    cout << "LU_factor() failed." << endl;
    return ans;
  }

  if (LU_solve(T, index, ans) != 0)  {
    cout << "LU_Solve() failed." << endl;
    return ans;
  }
  return ans;
} 

DMatrix solve(const DMatrix &a) {
  DMatrix b = ident(a.dim(1));
  return solve(a, b);
}

// transp(A) * inv(B) * C
template <class Matrix, class T>
Fortran_Matrix<T> matmult(
    const Transpose_View<Matrix> & A, 
    const Fortran_Matrix<T> &B)
{
    Subscript  M = A.num_rows();
    Subscript  N = A.num_cols();

    assert(B.num_rows() == N);
    Subscript L = B.num_cols();

    Fortran_Matrix<T> x(M,L);
    Subscript i, j, k;
    T tmp = 0;

    for (i=1; i<=M; i++) {
      for (j=1; j<=L; j++) {
        tmp = 0;
	for (k = 1; k <= N; k++) tmp += A(i,k) * B(k,j);
	x(i,j) = tmp;
      }
    }

    return x;
}

template <class Matrix, class T>
inline Fortran_Matrix<T> operator*(const Transpose_View<Matrix> & A, const Fortran_Matrix<T> &B)
{
    return matmult(A,B);
}



DMatrix AtBiC(const DMatrix &A, const DMatrix &B, const DMatrix &C) {
  DMatrix BiC = solve(B, C);
  return Transpose_view(A) * BiC;
}

DVector AtBiC(const DMatrix &A, const DMatrix &B, const DVector &C) {
  DVector BiC = solve(B, C);
  return Transpose_view(A) * BiC;
}

//crossprod
template <class T>
Fortran_Matrix<T> crossprod(const Vector<T> &v) {
  int n = v.size();
  Fortran_Matrix<T> ans(n,n);
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      ans(i,j) = v(i) * v(j);
  return ans;
}

template <class T>
Fortran_Matrix<T> crossprod(const Vector<T> &v1, const Vector<T> &v2) {
  int m = v1.size(), n = v2.size();
  Fortran_Matrix<T> ans(m,n);
  for (int i = 1;  i <= m; i++)
    for (int j = 1; j <= n; j++)
      ans(i,j) = v1(i) * v2(j);
  return ans;
}
