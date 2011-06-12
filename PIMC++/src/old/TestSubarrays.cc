#include "Common/Blitz.h"

main()
{
  Array<double,1> A(10);

  A = tensor::i;
  cerr << A << endl;
  Array<double,1> B;
  B.reference (A(Range(0,8,2)));
  cerr << B << endl;

  B = tensor::i;
  cerr << A << endl;


  Array<double,2> C(10,10);
  C = 10*tensor::i + tensor::j;
  Array<double,2> D;
  D.reference(C(Range::all(),Range(3,5)));

  cerr << "D = " << endl << D << endl;

}
