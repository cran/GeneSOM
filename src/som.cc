#include <string>

#include "SomParam.h"

Subscript find_winner(DMatrix &data, Subscript obs,
		      DMatrix &code) {
  Subscript winner = 1;
  Index1D C(1, data.num_cols());
  Index1D R(1, 1);
  Index1D Obs(obs, obs);
  double min_dist = norm2(data(Obs, C) - code(R, C));
  
  for (Subscript i = 2; i <= code.num_rows(); i++) {
    Index1D R(i, i);
    double dd = norm2(data(Obs, C) - code(R, C));
    if (dd < min_dist) {
      winner = i;
      min_dist = dd;
    }
  }
  return winner;
}

int update(DMatrix &code, DMatrix &data,
	   Subscript obs, double alpha, DVector &neigh) {
  Index1D Obs(obs, obs), Cols(1, data.num_cols());
  DMatrix dlt(code.num_rows(), code.num_cols());
  for (int i = 1; i <= code.num_rows(); i++) {
    Index1D I(i, i);
    dlt(I, Cols) = alpha * neigh(i) * 
      (data(Obs, Cols) - code(I, Cols));
  }
  code = code + dlt;
  return 0;
}

DMatrix GenCord(int xdim, int ydim) {
  DMatrix cord(xdim * ydim, 2);
  for (int i = 1; i <= ydim; i++) {
    for (int j = 1; j <= xdim; j++) {
      int r = xdim * (i - 1) + j;
      cord(r, 1) = j - 1;
      cord(r, 2) = i - 1;
    }
  }
  return cord;
}

void visual(DMatrix &data, DMatrix &code,
	    DMatrix &cord, DMatrix &vis) {
  Index1D Cols(1, data.num_cols());
  for (Subscript i = 1; i <= data.num_rows(); i++) {
    Subscript winner = find_winner(data, i, code);
    Index1D I(i,i), J(1,2), Win(winner, winner);
    vis(I, J) = cord(Win, J);
    double dd = norm2(data(I,Cols) - code(Win, Cols));
    vis(i, 3) = sqrt(dd);
  }
}

double qerror(DMatrix &data, DMatrix &code, 
	      DMatrix &cord, DMatrix &vis, SomParam &p) {
  int size = p.mapsize();
  double qerr = .0;
  Index1D Cols(1, data.num_cols());

  for (Subscript i = 1; i <= data.num_rows(); i++) {
    Subscript winner = cord2index(vis(i, 1), vis(i, 2), p.xdim());
    DVector nei = p.neigh(cord, winner, p.qerror_radius());

    Index1D I(i,i);
    for (int j = 1; j <= size; j++) {
      if (nei(j) == 0.0) continue;
      Index1D J(j,j);
      double dd = norm2(data(I,Cols) - code(J, Cols));
      qerr += nei(j) * dd;      
    }
  }
  return qerr / data.num_rows();
}

void som_train(DMatrix &data, DMatrix &code, DMatrix &cord, 
	       DMatrix &vis, SomParam &p) {
  //int size = p.mapsize();
  for (Subscript i = 1; i <= p.rlen(); i++) {
    Subscript obs = (i - 1) % data.num_rows() + 1;
    Subscript winner = find_winner(data, obs, code);

    double alp = p.alpha(i);
    double rad = p.radius(i);

    DVector nei = p.neigh(cord, winner, rad);

    update(code, data, obs, alp, nei);
  }
}

void som_top(DMatrix &data, DMatrix &code, DMatrix &vis, 
	     SomParam &p1, SomParam &p2, 
	     double *qerr) {
  // p1 and p2 only differs in alpha0 and radius0, maybe rlen;
  DMatrix cord = GenCord(p1.xdim(), p1.ydim());
  som_train(data, code, cord, vis, p1);
  som_train(data, code, cord, vis, p2);
  visual(data, code, cord, vis);
  *qerr = qerror(data, code, cord, vis, p1);
}


extern "C"{
#include <R.h>
#include <Rdefines.h>
}

DMatrix asDMatrix(SEXP a) {
  double *x;
  x = NUMERIC_POINTER(a);
  int *dims = INTEGER_POINTER(AS_INTEGER(GET_DIM(a)));
  DMatrix ans(dims[0], dims[1], x);
  return ans;
}

SEXP asSEXP(const DMatrix &a) {
  int size = a.num_cols() * a.num_rows();

  SEXP val;
  PROTECT(val = NEW_NUMERIC(size));
  double *p = NUMERIC_POINTER(val);
  const double *q = a.begin();
  for (int i = 0; i < size; i++) p[i] = q[i];
  SET_CLASS(val, ScalarString(mkChar("matrix")));

  SEXP dim;
  PROTECT(dim = NEW_INTEGER(2));
  INTEGER(dim)[0] = a.num_rows(); INTEGER(dim)[1] = a.num_cols();
  SET_DIM(val, dim);

  UNPROTECT(2);
  return val;
}

SEXP getListElement(SEXP list, char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

SomParam asSomParam(SEXP p) {
  int AlphaType, NeighType, TopolType, xdim, ydim, rlen;
  double alpha0, radius0, err_radius, inv_alp_c;
  AlphaType = INTEGER(VECTOR_ELT(p, 0))[0];
  NeighType = INTEGER(VECTOR_ELT(p, 1))[0];
  TopolType = INTEGER(VECTOR_ELT(p, 2))[0];
  alpha0 = REAL(VECTOR_ELT(p, 3))[0];
  radius0 = REAL(VECTOR_ELT(p, 4))[0];
  rlen = (int) REAL(VECTOR_ELT(p, 5))[0];
  err_radius = REAL(VECTOR_ELT(p, 6))[0];
  xdim = (int) REAL(VECTOR_ELT(p, 7))[0];
  ydim = (int) REAL(VECTOR_ELT(p, 8))[0];
  inv_alp_c = REAL(VECTOR_ELT(p, 9))[0];
  SomParam par(AlphaType, NeighType, TopolType, 
	       alpha0, radius0, rlen,
	       err_radius, xdim, ydim, inv_alp_c);
  return par;
}

extern "C"{
  SEXP som(SEXP data, SEXP code, SEXP p) {
    DMatrix Data = asDMatrix(data), Code = asDMatrix(code);
    SomParam P = asSomParam(p);
    DMatrix Cord = GenCord(P.xdim(), P.ydim());
    DMatrix Vis(Data.num_rows(), 3);
    som_train(Data, Code, Cord, Vis, P);
    visual(Data, Code, Cord, Vis);
    double qerr = qerror(Data, Code, Cord, Vis, P);
    SEXP ans, qe, names;
    PROTECT(ans = NEW_LIST(3));
    SET_VECTOR_ELT(ans, 0, asSEXP(Code));
    SET_VECTOR_ELT(ans, 1, asSEXP(Vis));
    PROTECT(qe = NEW_NUMERIC(1));
    NUMERIC_POINTER(qe)[0] = qerr;
    SET_VECTOR_ELT(ans, 2, qe);    

    PROTECT(names = NEW_STRING(3));
    SET_STRING_ELT(names, 0, mkChar("code"));
    SET_STRING_ELT(names, 1, mkChar("visual"));
    SET_STRING_ELT(names, 2, mkChar("qerror"));
    
    SET_NAMES(ans, names);
    UNPROTECT(3);
    return ans;
  }

  SEXP som_bat(SEXP data, SEXP code, SEXP param1, SEXP param2) {
    DMatrix Data = asDMatrix(data), Code = asDMatrix(code);
    //cout << "code = " << Code;
    SomParam p1 = asSomParam(param1);
    SomParam p2 = asSomParam(param2);
    double qerr = 0.0;
    DMatrix Vis(Data.num_rows(), 3);
    som_top(Data, Code, Vis, p1, p2, &qerr);

    SEXP ans, qe, names;
    PROTECT(ans = NEW_LIST(3));
    SET_VECTOR_ELT(ans, 0, asSEXP(Code));
    SET_VECTOR_ELT(ans, 1, asSEXP(Vis));
    PROTECT(qe = NEW_NUMERIC(1));
    NUMERIC_POINTER(qe)[0] = qerr;
    SET_VECTOR_ELT(ans, 2, qe);    

    PROTECT(names = NEW_STRING(3));
    SET_STRING_ELT(names, 0, mkChar("code"));
    SET_STRING_ELT(names, 1, mkChar("visual"));
    SET_STRING_ELT(names, 2, mkChar("qerror"));
    
    SET_NAMES(ans, names);
    UNPROTECT(3);
    return ans;
  }
  /*
  SEXP qerror_rap(SEXP data, SEXP code, SEXP vis, SEXP p) {
    DMatrix Data = asDMatrix(data), Code = asDMatrix(code), Vis = asDMatrix(vis);
    SomParam P = asSomParam(p);
    DMatrix Cord = GenCord(p.xdim(), p.ydim());
    double qerr = qerror(Data, Code, Cord, Vis);
    SEXP ans;
    PROTECT(ans = NEW_NUMERIC(1));
    NUMERIC_POINTER(ans)[0] = qerr;
    UNPROTECT(1);
    return ans;
  }
  */
}
