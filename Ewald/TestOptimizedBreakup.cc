#include "OptimizedBreakup.h"

#include <cstdio>

void WriteBasis(BasisClass &basis)
{
  FILE *fout = fopen ("LPQHI.dat", "w");
  for (double r=0.0; r<=basis.Get_rc(); r+=0.01) {
    for (int n=0; n<basis.NumElements(); n++)
      fprintf (fout, "%1.12e ", basis.h(n, r));
    fprintf (fout, "\n");
  }
  fclose (fout);
}

void  TestLPQHI()
{
  LPQHI_BasisClass basis;

  basis.SetNumKnots(6);
  basis.Set_rc(5.0);

  WriteBasis (basis);
}


main()
{
  TestLPQHI();
}
