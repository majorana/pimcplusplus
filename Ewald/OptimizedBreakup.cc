#include "OptimizedBreakup.h"

/////////////////////////////////////////
/// LPQHI_BasisClass Member Functions ///
/////////////////////////////////////////

void LPQHI_BasisClass::SetNumKnots(int n)
{
  assert (n > 1);
  NumKnots = n;
  if (r_c != 0.0) {
    delta = r_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
}

int LPQHI_BasisClass::NumElements()
{
  return 3*NumKnots;
}

void LPQHI_BasisClass::Set_rc(double rc)
{
  r_c = rc;
  if (NumKnots != 0) {
    delta = r_c / (NumKnots - 1);
    deltaInv = 1.0/delta;
  }
  cerr << "delta = " << delta << endl;
}


double LPQHI_BasisClass::h(int n, double r)
{
  int i=n/3;
  int alpha = n-3*i;
  double ra = delta*(i-1);
  double rb = delta*i;
  double rc = delta*(i+1);
  if ((r > ra) && (r <= rb)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += (S(alpha,j) * prod);
      prod *= ((rb - r) * deltaInv);
    }
    for (int j=0; j<alpha; j++)
      sum *= (-delta);
    return (sum);
  }
  else if ((r > rb) && (r <= rc)) {
    double sum = 0.0;
    double prod = 1.0;
    for (int j=0; j<=5; j++) {
      sum += S(alpha,j) * prod;
      prod *= ((r-rb) * deltaInv);
    }
    for (int j=0; j<alpha; j++)
      sum *= delta;
    return sum;
  }
  return 0.0;
};


double LPQHI_BasisClass::c(int n, double k)
{
  int i=n/3;
  int alpha = n-3*i;

};
