#include <math.h>
#include <iostream.h>

#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H
#include "blas++.h"

#define INV_ALPHA_CONST 100.0

void rect2hexa(double x, double y, double *u, double *v) {
  if (((int)(y) % 2) == 0) {
    *u = x;
    *v = sqrt(3) / 2 * y;
  }
  else {
    *u = x + 1.0/2;
    *v = sqrt(3) / 2 * y;
  }
}

double norm2(LaVectorDouble x) {
  double ans = 0.0;
  for (int i = 0; i < x.size(); i++) 
    ans += x(i)*x(i);
  ans = sqrt(ans);
  return ans;
}

double norm2(LaGenMatDouble x) {
  double ans = 0.0;
  for (int i = 0; i < x.size(0); i++) 
    for (int j = 0; j < x.size(1); j++)
      ans += x(i, j)*x(i, j);
  ans = sqrt(ans);
  return ans;
}

LaVectorDouble rect2hexa(LaVectorDouble &x) {
  double u[2];
  rect2hexa(x(0), x(1), &u[0], &u[1]);
  LaVectorDouble y(u, 2);
  return y;
}

LaGenMatDouble rect2hexa(LaGenMatDouble &x) {
  double u[2];
  rect2hexa(x(0, 0), x(0, 1), &u[0], &u[1]);
  LaGenMatDouble y(u, 1, 2);
  return y;
}


double lin_radius(double radius_0, int t, int rlen) {
  return 1.0 + (radius_0 - 1.0) * (double) (rlen - t) / (double) rlen;
}


int find_winner(LaGenMatDouble &data, int size, int obs, 
		LaGenMatDouble &code) {
  int winner = 0;
  double min_dist = norm2(data(LaIndex(obs), LaIndex()) - 
			       code(LaIndex(0), LaIndex()));
  for (int i = 1; i < size; i++) {
    double dd = norm2(data(LaIndex(obs), LaIndex()) - 
		    code(LaIndex(i), LaIndex()));
    if (dd < min_dist) {
      winner = i;
      min_dist = dd;
    }
  }
  return winner;
} 

LaGenMatDouble SMult(double v, LaGenMatDouble &m) {
  LaGenMatDouble tmp = m;
  for (int i = 0; i < m.size(0); i++)
    for (int j = 0; j < m.size(1); j++) {
      tmp(i, j) = v * m(i, j);
    }
  return tmp;
}

LaGenMatDouble Smult(LaVectorDouble &v, LaGenMatDouble &m) {
  LaGenMatDouble tmp = m;
  for (int i = 0; i < m.size(0); i++)
    for (int j = 0; j < m.size(1); j++)
      tmp(i, j) = v(i) * m(i, j);
  return tmp;
}

int update(LaGenMatDouble &code, int size, 
	   LaGenMatDouble &data, int obs, 
	   double alpha, LaVectorDouble neigh) {
  for (int i = 0; i < size; i++) {
    LaGenMatDouble tmp = data(LaIndex(obs), LaIndex()) - code(LaIndex(i), LaIndex());
    tmp = SMult(alpha * neigh(i),  tmp);
    tmp = code(LaIndex(i), LaIndex()) + tmp;
    code(LaIndex(i), LaIndex()) = tmp;
    //code(LaIndex(i), LaIndex()) = code(LaIndex(i), LaIndex()) + 
    //  alpha * neigh(i) * tmp;
    //      (data(LaIndex(obs), LaIndex()) - code(LaIndex(i), LaIndex()));
  }
  return 0;
}

LaGenMatDouble genCord(int xdim, int ydim) {
  LaGenMatDouble cord(xdim * ydim, 2);
  for (int i = 0; i < xdim; i++) { 
    for (int j = 0; j < ydim; j++) {
      cord(LaIndex(ydim * i + j), LaIndex(0)) = i;
      cord(LaIndex(ydim * i + j), LaIndex(1)) = j;
    }
  }
  return cord;
}


class SomParam {
public:
  typedef double SomParam::Alpha(double alpha_0, int t, int rlen);
  typedef double SomParam::Dist(LaGenMatDouble v1, LaGenMatDouble v2);
  typedef LaVectorDouble SomParam::Neigh(LaGenMatDouble &cord,
			       int size, int winner, 
			       double radius, 
			       double (*Dist)(LaGenMatDouble, LaGenMatDouble));

private:
  static double rect_dist(LaGenMatDouble v1, LaGenMatDouble v2) {
    return norm2(v1 - v2);
  };

  static double hexa_dist(LaGenMatDouble v1, LaGenMatDouble v2) {
    LaGenMatDouble u1 = rect2hexa(v1);
    LaGenMatDouble u2 = rect2hexa(v2);
    return rect_dist(u1, u2);
  };
  
  static double lin_alpha(double alpha_0, int t, int rlen) {
    return alpha_0 * (1.0 - (double) t / rlen);
  };
  
  static double inv_alpha(double alpha_0, int t, int rlen) {
    double C = rlen / INV_ALPHA_CONST; 
    return alpha_0 * C / (C + t);
  };
  
  static LaVectorDouble bubble_neigh(LaGenMatDouble &cord,
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
  
  static LaVectorDouble gaussian_neigh(LaGenMatDouble &cord,
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

  Alpha *_alpha;
  Dist *_dist;
  Neigh *_neigh;
  //static Alpha *Alpha_list[2] = {lin_alpha, inv_alpha};
  //static Dist *Dist_list[2] = {rect_dist, hexa_dist};
  //static Neigh *Neigh_list[2] = {bubble_neigh, gaussian_neigh};
public:

  SomParam(int alphaType, int neighType, int topol) {
    Alpha *Alpha_list[2] = {lin_alpha, inv_alpha};
    Dist *Dist_list[2] = {rect_dist, hexa_dist};
    Neigh *Neigh_list[2] = {bubble_neigh, gaussian_neigh};
    _alpha = Alpha_list[alphaType - 1];
    _dist = Dist_list[topol - 1];
    _neigh = Neigh_list[neighType - 1];
    //if (alphaType == 1) _alpha = lin_alpha;
    //else _alpha = inv_alpha;
    //if (topol == 1) _dist = rect_dist;
    //else _dist = hexa_dist;
    //if (neighType == 1) _neigh = bubble_neigh;
    //else _neigh = gaussian_neigh;
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
  //double (*Dist)(LaGenMatDouble, LaGenMatDouble)) {
    return _neigh(cord, size, winner, radius, _dist);
  } 

  ~SomParam() {}
};


void visual(LaGenMatDouble &data, LaGenMatDouble &code, 
	    LaGenMatDouble &cord, int size, 
	    LaGenMatDouble &vis) {
  for (int i = 0; i < data.size(0); i++) {
    int winner = find_winner(data, size, i, code);
    vis(LaIndex(i), LaIndex(0, 1)) = cord(LaIndex(winner), LaIndex());
    double dd = norm2(data(LaIndex(i), LaIndex()) - code(LaIndex(winner), LaIndex()));
    vis(i, 2) = dd * dd;
  }
}

int cord2index(double x, double y, int xdim) {
  return (int) (x + y * xdim);
}

double qerror(LaGenMatDouble &data, LaGenMatDouble &code, 
	      LaGenMatDouble &cord, int xdim, int ydim, 
	      LaGenMatDouble &vis, SomParam &p,
	      double radius) {
  int size = xdim * ydim;
  double qerr = 0;
  for (int i = 0; i < data.size(0); i++) {
      int winner = cord2index(vis(i, 0), vis(i, 1), xdim);
      LaVectorDouble nei = p.neigh(cord, size, winner, radius);
      for (int j = 0; j < size; j++) {
	if (nei(j) == 0) continue;
	double dd = norm2(data(LaIndex(i), LaIndex()) - code(LaIndex(j), LaIndex()));
	qerr += nei(j) * dd * dd;      
    }
  }
  return qerr;
}

extern "C" {  
  void som(double *data, 
	   int *nrow, int *ncol,
	   double *code,
	   int *xdim, int *ydim,
	   double *alpha, int *alphaType,
	   int *neigh, int *topol,
	   double *radius, int *rlen,
	   double *vis) {
    LaGenMatDouble Data(data, *nrow, *ncol);

    int size = *xdim * *ydim;
    LaGenMatDouble Code(code, size, *ncol);
    //cout << "Initial Code book: " << endl << Code;

    //set up cord matrix
    LaGenMatDouble Cord = genCord(*xdim, *ydim);

    SomParam p(*alphaType, *neigh, *topol);

    double (*Radius)(double radius_0, int t, int rlen);
    Radius = lin_radius;
    
    for (int train = 0; train < 2; train++) {
      for (int i = 0; i < rlen[train]; i++) {
	int obs = i % *nrow;
	int winner = find_winner(Data, size, obs, Code);

	//get alpha, radius, neigh
	double alp = p.alpha(alpha[train], i, rlen[train]);
	double rad = Radius(radius[train], i, rlen[train]);
	
	LaVectorDouble nei = p.neigh(Cord, size, winner, rad);
	//the dist type is hidden in the p now as p._dist;

	update(Code, size, Data, obs, alp, nei);
      }
    }
    // visualization
    LaGenMatDouble Vis(vis, *nrow, 3);
    visual(Data, Code, Cord, size, Vis);
  }

  void find_qerror(double *data, 
		   int *nrow, int *ncol,
		   double *code,
		   int *xdim, int *ydim,
		   int *alphaType,
		   int *neigh, int *topol,
		   double *radius, 
		   double *vis, double *qerr) {
    LaGenMatDouble Data(data, *nrow, *ncol);

    int size = *xdim * *ydim;
    LaGenMatDouble Code(code, size, *ncol);
    //cout << "Initial Code book: " << endl << Code;
    
    LaGenMatDouble Vis(vis, size, 3);
    //cout << "Vis is :" << Vis;
    LaGenMatDouble Cord = genCord(*xdim, *ydim);
    //cout << "cord is " << endl << Cord;

    SomParam p(*alphaType, *neigh, *topol);

    *qerr = qerror(Data, Code, 
		   Cord, *xdim, *ydim, 
		   Vis, p,
		   *radius);
  }

}
