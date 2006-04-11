#include "OptimizedBreakup.h"

#include <cstdio>

void WriteBasis(BasisClass &basis)
{
  FILE *fout = fopen ("LPQHI.dat", "w");
  for (double r=0.0; r<=basis.Get_rc(); r+=0.001) {
    fprintf (fout, "%1.12e ", r);
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

void Test_kDeriv()
{
  LPQHI_BasisClass basis;
  basis.SetNumKnots(6);
  basis.Set_rc(5.0);
  double epsilon=1.0e-6;
  double k = 0.347;

  for (int i=0; i<5; i++)
    for (int n=0; n<5; n++) {
      complex<double> FD = 
	(basis.Dminus(i,k+epsilon,n)-basis.Dminus(i,k-epsilon,n))/(2.0*epsilon);
      complex<double> an = basis.dDminus_dk(i,k,n);
      fprintf (stderr, "%1.12e + %1.12ei      %1.12e + %1.12ei\n",
	       real(FD), imag(FD), real(an), imag(an));
    }
}

void TestBasisDeriv()
{
  LPQHI_BasisClass basis;
  basis.SetNumKnots(3);
  basis.Set_rc(2.0);
  double epsilon=1.0e-5;
  double k = 1.0;
  
  for (int i=0; i<basis.NumElements(); i++) {
    double FD = (basis.c(i,k+epsilon)-basis.c(i,k-epsilon))/(2.0*epsilon);
    double an = basis.dc_dk(i,k);
    fprintf (stderr, "%1.12e  %1.12e\n", FD, an);
  }
}







void TestCoulomb()
{
  LPQHI_BasisClass basis;
  TinyVector<double,3> box;
  box = 10.0, 10.0, 10.0;
  basis.SetBox (box);
  double Omega = box[0]*box[1]*box[2];
  basis.SetNumKnots(20);
  basis.Set_rc(5.0);

//   double k = 1;
//   for (int n=0; n<basis.NumElements(); n++) {
//     double c = basis.c(n,k);
//     double c_numerical = basis.c_numerical(n,k);
//     fprintf (stderr, "n = %d\n", n);
//     fprintf (stderr, "c           = %1.12e\n", c);
//     fprintf (stderr, "c_numerical = %1.12e\n", c_numerical);
//   }

  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (2.0, 25.0, 1000.0);
  double delta = basis.GetDelta();

  Array<double,1> Vk(breakup.kpoints.size());
  for (int i=0; i<breakup.kpoints.size(); i++) {
    double k = breakup.kpoints(i)[0];
    double k0 = 2.0*M_PI;
    Vk(i) = 4.0*M_PI/(Omega*(k*k)) * cos(k*basis.Get_rc());
  }
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1> adjust(N);
  adjust = true;
  t = 0.0;
  adjust(N-3) = false; t(N-3) = -0.2;
  adjust(N-2) = false; t(N-2) = 1.0/25.0*delta;
  adjust(N-1) = false; t(N-1) = -1/125.0*delta*delta;
  adjust(1)   = false; t(1) = 0.0;

  double chi2 = breakup.DoBreakup(Vk, t, adjust);
  cerr << "chi-squared = " << chi2 << endl;
  cerr << "t = " << t << endl;

  FILE *fout = fopen ("Vlong.dat", "w");
  
  for (double r=0.0; r<10.0; r+=0.001) {
    double v = 0.0;
    for (int n=0; n<t.rows(); n++)
      v += t(n)*basis.h(n,r);
    fprintf (fout, "%1.12e %1.12e\n", r, v);

  }
  fclose(fout);

}


main()
{
  TestBasisDeriv();
//   TestLPQHI();
//   TestCoulomb();
}
