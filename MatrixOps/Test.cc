#include "MatrixOps.h"

main()
{
  Array<double,2> A(3,3);
  A = 1.2,  -5.6,  7.3,
     -2.1,  -9.5,  7.1,
      4.9,   1.6, -6.6;

  Array<double,1> perm;
  double sign;
  cerr << "A = " << A << endl;
  LUdecomp (A, perm, sign);
  cerr << "LU = " << A << endl;

}
