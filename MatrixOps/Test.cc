#include "MatrixOps.h"

main()
{
  Array<double,2> A(3,3), LU(3,3), Ainv(3,3);
  A = 1.2,  -5.6,  7.3,
     -2.1,  -9.5,  7.1,
      4.9,   1.6, -6.6;
  Array<double,1> b(3), x(3);
  
  b = 1.0, 2.0, 3.0;

  LU = A;
  Array<int,1> perm;
  double sign;
  //cerr << "A = " << A << endl;
  LUdecomp (LU, perm, sign);
  x = b;
  LUsolve (LU, perm, x);
  //cerr << "LU = " << LU << endl;
  //cerr << "b = " << b << endl;
  //cerr << "x = " << x << endl;
  Ainv = Inverse(A);
  //cerr << "Ainv = " << Ainv << endl;
  
  Array<double,2> AinvA(3,3);
  AinvA = 0.0;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
	AinvA(i,j) += Ainv(i,k)*A(k,j);
  bool success = true;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      if ((i==j)) {
	if (fabs(AinvA(i,j) - 1.0) > 1e-14)
	  success = false;
      }
      else {
	if (fabs(AinvA(i,j)) > 1e-14)
	  success = false;
      }
  cerr << (success ? "Passed Matrix Inverse Test\n" : 
    "Failed Matrix Inverse Test\n");

  // SVD test
  Array<double,2> U(3,3), V(3,3);
  Array<double,1> S(3);
  Array<double,2> B(4,3);

  B = 1.2,  -5.6,  7.3,
     -2.1,  -9.5,  7.1,
      4.9,   1.6, -6.6,
      1.3,   1.8,  9.0;

  SVdecomp (B, U, S, V);
  cerr << "B = " << B << endl;
  cerr << "U = " << U << endl;
  cerr << "V = " << V << endl;
  cerr << "S = " << S << endl;

}
