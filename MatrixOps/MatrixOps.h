#include "../Blitz.h"

void LUdecomp (Array<double,2> A, Array<double,1> perm, 
	       double &sign);

void LUsolve (Array<double,2> &LU, Array<double,1> perm,
	      Array<double,1> &b);

