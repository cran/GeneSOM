void rect2hexa(double x, double y, double *u, double *v);

double norm2(LaVectorDouble x);

double norm2(LaGenMatDouble x);

LaVectorDouble rect2hexa(LaVectorDouble &x);

LaGenMatDouble rect2hexa(LaGenMatDouble &x);

double lin_radius(double radius_0, int t, int rlen);

int find_winner(LaGenMatDouble &data, int size, int obs, 
		LaGenMatDouble &code);

LaGenMatDouble SMult(double v, LaGenMatDouble &m);

LaGenMatDouble Smult(LaVectorDouble &v, LaGenMatDouble &m);

int update(LaGenMatDouble &code, int size, 
	   LaGenMatDouble &data, int obs, 
	   double alpha, LaVectorDouble neigh);

LaGenMatDouble genCord(int xdim, int ydim);

void visual(LaGenMatDouble &data, LaGenMatDouble &code, 
	    LaGenMatDouble &cord, int size, 
	    LaGenMatDouble &vis);

int cord2index(double x, double y, int xdim);
