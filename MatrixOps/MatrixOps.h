#include "../Blitz.h"

void LUdecomp (Array<double,2> &A, Array<int,1> &perm, 
	       double &sign);

void LUsolve (Array<double,2> &LU, Array<int,1> &perm,
	      Array<double,1> &b);

Array<double,2> Inverse (Array<double,2> &A);
