#include <math.h>
#include <iostream.h>

#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H
#include "blas++.h"

#define INV_ALPHA_CONST 100.0

#include "som.h"

class SomParam {
public:
  typedef double SomParam::Alpha(double alpha_0, int t, int rlen);
  typedef double SomParam::Dist(LaGenMatDouble v1, LaGenMatDouble v2);
  typedef LaVectorDouble SomParam::Neigh(LaGenMatDouble &cord,
			       int size, int winner, 
			       double radius, 
			       double (*Dist)(LaGenMatDouble, LaGenMatDouble));
protected:
  inline static double rect_dist(LaGenMatDouble v1, LaGenMatDouble v2);
  inline static double hexa_dist(LaGenMatDouble v1, LaGenMatDouble v2);
  inline static double lin_alpha(double alpha_0, int t, int rlen);
  inline static double inv_alpha(double alpha_0, int t, int rlen);
  inline static LaVectorDouble bubble_neigh(LaGenMatDouble &cord,
				     int size, int winner, 
				     double radius,
				     double (*dist)(LaGenMatDouble, 
						    LaGenMatDouble));
  inline static LaVectorDouble gaussian_neigh(LaGenMatDouble &cord,
				       int size, int winner, 
				       double radius,
				       double (*dist)(LaGenMatDouble, 
						      LaGenMatDouble));
  Alpha *_alpha;
  Dist *_dist;
  Neigh *_neigh;

public:

  SomParam(int alphaType, int neighType, int topol) {
    Alpha *Alpha_list[2] = {lin_alpha, inv_alpha};
    Dist *Dist_list[2] = {rect_dist, hexa_dist};
    Neigh *Neigh_list[2] = {bubble_neigh, gaussian_neigh};
    _alpha = Alpha_list[alphaType - 1];
    _dist = Dist_list[topol - 1];
    _neigh = Neigh_list[neighType - 1];
  }
  double alpha(double alpha_0, int t, int rlen) {
    return _alpha(alpha_0, t, rlen);
  }
  double dist(LaGenMatDouble v1, LaGenMatDouble v2) {
    return _dist(v1, v2);
  }
  LaVectorDouble neigh(LaGenMatDouble &cord,
		       int size, int winner, 
		       double radius) {
    return _neigh(cord, size, winner, radius, _dist);
  } 

  ~SomParam() {}
};

double SomParam::rect_dist(LaGenMatDouble v1, LaGenMatDouble v2) {
  return norm2(v1 - v2);
};

double SomParam::hexa_dist(LaGenMatDouble v1, LaGenMatDouble v2) {
    LaGenMatDouble u1 = rect2hexa(v1);
    LaGenMatDouble u2 = rect2hexa(v2);
    return rect_dist(u1, u2);
};
  
double SomParam::lin_alpha(double alpha_0, int t, int rlen) {
    return alpha_0 * (1.0 - (double) t / rlen);
};
  
double SomParam::inv_alpha(double alpha_0, int t, int rlen) {
  double C = rlen / INV_ALPHA_CONST; 
  return alpha_0 * C / (C + t);
};
  
LaVectorDouble SomParam::bubble_neigh(LaGenMatDouble &cord,
				     int size, int winner, 
				     double radius,
				     double (*dist)(LaGenMatDouble, 
						    LaGenMatDouble)) {
    int i;
    double dd;
    LaVectorDouble neigh(size);
    for (i = 0; i < size; i++) {
      dd = dist(cord(LaIndex(i), LaIndex()),
		      cord(LaIndex(winner), LaIndex()));
      neigh(i) = (dd <= radius) ? 1.0 : 0.0;
    }
    return neigh;
};
  
LaVectorDouble SomParam::gaussian_neigh(LaGenMatDouble &cord,
					       int size, int winner, 
					       double radius,
					       double (*dist)(LaGenMatDouble, 
							      LaGenMatDouble)) {
  int i;
  double dd;
  LaVectorDouble neigh(size);
    for (i = 0; i < size; i++) {
      dd = dist(cord(LaIndex(i), LaIndex()),
		cord(LaIndex(winner), LaIndex()));
      neigh(i) = exp( - dd * dd / 2.0 / radius / radius);
    }
    return neigh;
};
