#include <math.h>
#include <iostream.h>

#include "lafnames.h"
#include LA_VECTOR_DOUBLE_H
#include "blas++.h"

#include "SomParam.h"

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
