#include "OptimizedBreakup.h"
#include "../MatrixOps/MatrixOps.h"
#include "../Integration/GKIntegration.h"

class cIntegrand
{
private:
  BasisClass &Basis;
  int n;
  double k;
public:
  double operator()(double r)
  {
    return 4.0*M_PI/(k*Basis.Omega) *Basis.h(n,r)*r*sin(k*r);
  }
  cIntegrand(BasisClass &basis, int n_, double k_) : Basis(basis), n(n_), k(k_)
  { /* do nothing */ }
};

///////////////////////////////////
/// BasisClass Member Functions ///
///////////////////////////////////
double BasisClass::c_numerical(int n, double k)
{
  int i = n/3;
  int numKnots = NumElements()/3;
  double ra = 0.0;
  double rb = r_c;

  cIntegrand integrand(*this, n, k);
  GKIntegration<cIntegrand,GK31> integrator(integrand);
  integrator.SetRelativeErrorMode();
  return integrator.Integrate (ra, rb, 1.0e-12, 1.0e-10, false);
}



///////////////////////////////////////////
/// OptimizedBreakup Memember Functions ///
///////////////////////////////////////////
void OptimizedBreakup::SetkVecs(double kc, double kMax)
{
  int numk = 0;
  TinyVector<double,3> b;
  b[0] = 2.0*M_PI/Basis.GetBox()[0];
  b[1] = 2.0*M_PI/Basis.GetBox()[1];
  b[2] = 2.0*M_PI/Basis.GetBox()[2];
  TinyVector<int,1> maxIndex;
  maxIndex[0] = (int)ceil(kMax/b[0]);
  maxIndex[1] = (int)ceil(kMax/b[1]);
  maxIndex[2] = (int)ceil(kMax/b[2]);

  TinyVector<double,3> k;
  for (int ix=-maxIndex[0]; ix<=maxIndex[0]; ix++) {
    k[0] = ix*b[0];
    for (int iy=-maxIndex[1]; iy<=maxIndex[1]; iy++) {
      k[1] = iy*b[1];
      for (int iz=-maxIndex[2]; iz<=maxIndex[2]; iz++) {
	k[2] = iz*b[2];
	double k2 = dot(k,k);
	if ((k2 > (kc*kc)) && (k2 <= (kMax*kMax)))
	  numk++;
      }
    }
  }
  
  kVecs.resize(numk);
  numk = 0;
  for (int ix=-maxIndex[0]; ix<=maxIndex[0]; ix++) {
    k[0] = ix*b[0];
    for (int iy=-maxIndex[1]; iy<=maxIndex[1]; iy++) {
      k[1] = iy*b[1];
      for (int iz=-maxIndex[2]; iz<=maxIndex[2]; iz++) {
	k[2] = iz*b[2];
	double k2 = dot(k,k);
	if ((k2 > (kc*kc)) && (k2 <= (kMax*kMax))) {
	  kVecs(numk) = k;
	  numk++;
	}
      }
    }
  }
}

void OptimizedBreakup::DoBreakup(const Array<double,1> &Vk, Array<double,1> &t,
				 const Array<bool,1> &adjust)
{
  const double tolerance = 1.0e-12;
  assert(t.rows()==adjust.rows());
  assert(t.rows()==Basis.NumElements());
  Array<double,2> A;
  Array<double,1> b;
  Array<double,2> cnk;

  int numElem = t.rows();
  /// HACK:  for simplicity, we ignore adjust for now
  //  for (int i=0; i<adjust.rows(); i++)
  //   if (!adjust(i))
  //    numElem--;
  A.resize(numElem, numElem);
  b.resize(numElem);
  cnk.resize(numElem,kVecs.rows());

  // Fill in cnk.
  for (int n=0; n<t.rows(); n++) {
    for (int ki=0; ki<kVecs.rows(); ki++) {
      double k = sqrt(dot(kVecs(ki),kVecs(ki)));
      cnk(n,ki) = Basis.c_numerical(n,k);
    }
  }

  // Now, fill in A and b
  A = 0.0;
  b = 0.0;
  for (int l=0; l<numElem; l++) {
    for (int ki=0; ki<kVecs.rows(); ki++) {
      b(l) += Vk(ki) * cnk(l, ki);
      for (int n=0; n<numElem; n++) 
	A(l,n) += cnk(l,ki)*cnk(n,ki);
    }
  }

  //  cerr << "A = " << A << endl;
  //cerr << "b = " << b << endl;

  // Now do SVD decomposition:
  Array<double,2> U(numElem, numElem), V(numElem, numElem);
  Array<double,1> S(numElem), Sinv(numElem);
  SVdecomp(A, U, S, V);
  
  // Zero out near-singular values
  double Smax=S(0);
  for (int i=1; i<S.size(); i++)
    Smax = max (S(i),Smax);
  for (int i=0; i<S.size(); i++)
    Sinv(i) = (S(i) < (tolerance*Smax)) ? 0.0 : (1.0/S(i));

  t = 0.0;
  // Compute t_n, checking singular values
//   for (int i=0; i<numElem; i++) {
//     double coef = 0.0;
//     for (int j=0; j<numElem; j++)
//       coef += U(j,i) * b(j);
//     coef *= Sinv(i);
//     for (int j=0; j<numElem; j++)
//       t(j) += coef * V(j,i);
//   }
  for (int k=0; k<numElem; k++)
    for (int i=0; i<numElem; i++) {
      double coef = 0.0;
      for (int j=0; j<numElem; j++)
	coef += U(j,i) * b(j);
      coef *= Sinv(i);
      t(k) += coef * V(k,i);
  }

      
  

}

void OptimizedBreakup::DoBreakup(const Array<double,1> &Vk, Array<double,1> &t)
{
  Array<bool,1> adjust(t.rows());
  adjust = true;
  DoBreakup (Vk, t, adjust);
}

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


double LPQHI_BasisClass::c(int m, double k)
{
  int i=m/3;
  int alpha = m-3*i;
  
  double sum = 0.0;
  if (i == 0) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n));
    }
  else if (i == (NumKnots-1)) 
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dminus(i,k,n)*sign);
    }
  else
    for (int n=0; n<=5; n++) {
      double sign = ((alpha+n)&1) ? -1.0 : 1.0;
      sum += S(alpha, n) * (Dplus(i,k,n) + Dminus(i,k,n)*sign);
    }
  for (int j=0; j<alpha; j++)
    sum *= delta;

  return (sum);
};
