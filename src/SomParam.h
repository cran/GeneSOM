#include <math.h>
#include "sub.h"

//#define INV_ALPHA_CONST 100.0

double rect_dist(const DMatrix &v1, const DMatrix &v2) {
  return norm2(v1 - v2);
}

double rect_dist(const Region2D<DMatrix> v1, const Region2D<DMatrix> v2) {
  return norm2(v1 - v2);
}

double rect_dist(const_Region2D<DMatrix> v1, const_Region2D<DMatrix> v2) {
  return norm2(v1 - v2);
}

double hexa_dist(const DMatrix &v1, const DMatrix &v2) {
  DVector u1 = hexa2rect(v1);
  DVector u2 = hexa2rect(v2);
  return norm2(u1 - u2);
}

double hexa_dist(const Region2D<DMatrix> v1, const Region2D<DMatrix> v2) {
  DVector u1 = hexa2rect(v1);
  DVector u2 = hexa2rect(v2);
  return norm2(u1 - u2);
}

double hexa_dist(const_Region2D<DMatrix> v1, const_Region2D<DMatrix> v2) {
  DVector u1 = hexa2rect(v1);
  DVector u2 = hexa2rect(v2);
  return norm2(u1 - u2);
}

double lin_alpha(double alpha_0, int t, Subscript rlen, double) {
  return alpha_0 * (1.0 - (double) t / rlen);
}
  
double inv_alpha(double alpha_0, int t, Subscript rlen, double C) {
  //double C = rlen / INV_ALPHA_CONST; 
  return alpha_0 * C / (C + t);
}

double lin_radius(double radius_0, int t, int rlen) {
  return 1.0 + (radius_0 - 1.0) * (double) (rlen - t) / (double) rlen;
}

DVector bubble_neigh(const DMatrix &cord, Subscript winner, double radius,
		     double (*dist)(const_Region2D<DMatrix>, 
				    const_Region2D<DMatrix>)) {
  Subscript m = cord.num_rows(), n = cord.num_cols();
  DVector neigh(m);
  Index1D J(1,n), Win(winner,winner);
  for (Subscript i = 1; i <= m; i++) {
    Index1D I(i,i);
    double dd = dist(cord(I,J), cord(Win,J));
    neigh(i) = (dd <= radius) ? 1.0 : 0.0;
  }
  return neigh;
}

DVector gaussian_neigh(const DMatrix &cord, Subscript winner, double radius,
		       double (*dist)(const_Region2D<DMatrix>, 
				      const_Region2D<DMatrix>)) {
  Subscript m = cord.num_rows(), n = cord.num_cols();
  DVector neigh(m);
  Index1D J(1,n), Win(winner,winner);
  for (Subscript i = 1; i <= m; i++) {
    Index1D I(i,i);
    double dd = dist(cord(I,J), cord(Win,J));
    neigh(i) = exp(- dd / 2.0 / radius / radius);
  }
  return neigh;
}

#define LIN_ALPHA 1
#define INV_ALPHA 2

#define RECT 1
#define HEXA 2

#define BUBBLE 1
#define GAUSSIAN 2

class SomParam{
public:
  typedef double Alpha(double, int, int, double);
  typedef double Radius(double, int, int);
  typedef double Dist(const_Region2D<DMatrix>, const_Region2D<DMatrix>);
  typedef DVector Neigh(const DMatrix &, Subscript, double, Dist *);
protected:
  Alpha *_alpha;
  Radius *_radius;
  Dist *_dist;
  Neigh *_neigh;
  int _xdim;
  int _ydim;
  double _alpha0;
  Subscript _rlen;
  double _radius0;
  double _qerror_radius;
  double _inv_alp_c;

public: 
  SomParam(int AlphaType, int NeighType, int TopolType,
	   double alpha0, double radius0, Subscript rlen, 
	   double qerror_radius, int xdim, int ydim, double inv_alp_c) {
    if (AlphaType == LIN_ALPHA) _alpha = lin_alpha;
    else _alpha = inv_alpha;
    if (TopolType == RECT) _dist = rect_dist;
    else _dist = hexa_dist;
    if (NeighType == BUBBLE)  _neigh = bubble_neigh; 
    else _neigh = gaussian_neigh;

    _radius = lin_radius; 
    
    _xdim = xdim;
    _ydim = ydim;
    _rlen = rlen;
    _alpha0 = alpha0;
    _radius0 = radius0;
    _qerror_radius = qerror_radius;
    _inv_alp_c = inv_alp_c;
    
  }

  ~SomParam() {} 
  double alpha(int t) {
    return _alpha(_alpha0, t, _rlen, _inv_alp_c);
  }
  double dist(const_Region2D<DMatrix> v1, const_Region2D<DMatrix> v2) {
    return _dist(v1, v2);
  }
  DVector neigh(const DMatrix &cord, Subscript winner, double radius) {
    return _neigh(cord, winner, radius, _dist);
  }
  double radius(int t) {
    return _radius(_radius0, t, _rlen);
  }
  double alpha0() {return _alpha0;}
  double radius0() {return _radius0;}
  int xdim() {return _xdim;}
  int ydim() {return _ydim;}
  int mapsize() {return _xdim * _ydim;}
  Subscript rlen() {return _rlen;}
  double qerror_radius() {return _qerror_radius;}
}; 

