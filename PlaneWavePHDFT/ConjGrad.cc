#include "ConjGrad.h"


void ConjGrad::Setup()
{
  int N = H.GVecs.size();
  c.resize (N);
  cnext.resize (N);
  Hcnext.resize (N);
  Hc.resize (N);
  Phi.resize (N);
  LastPhi.resize (N);
  Phip.resize(N);
  Phipp.resize(N);
  Xi.resize(N);
  Xip.resize(N);
  Eta.resize(N);
  Etap.resize(N);
  c = 1.0;
  Normalize (c);
  IsSetup = true;
}

void ConjGrad::CalcPhiSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  cerr << "E = " << E0 << endl;
  Phip = E0*c - Hc;
  Normalize (Phip);  
}

void ConjGrad::Precondition()
{
  for (int i=0; i<c.size(); i++) {
    double x = 0.5*dot(H.GVecs(i), H.GVecs(i))/T;
    double num = 27.0 + 18.0*x +12.0*x*x + 8.0*x*x*x;
    double denom = num + 16.0*x*x*x*x;
    Eta(i) = (num/denom)* Xip(i);
  }
}

void ConjGrad::CalcPhiCG()
{
  Hc = 0.0;
  H.Kinetic.Apply (c, Hc);
  T = realconjdot (c, Hc);
  H.CoulombFFT.Apply (c, Hc);
  E0 = realconjdot (c, Hc);
  Xi = E0*c - Hc;
  /// Orthonalize to other bands here
  Xip = Xi;
  
  Precondition();
  //Eta = Xip;

  // Now, orthogonalize to psi
  Etap = Eta - conjdot (c, Eta)*c;
  complex<double> etaxi = conjdot(Etap, Xip);
  complex<double> gamma; 
  if (EtaXiLast != complex<double>(0.0, 0.0)) 
    gamma = etaxi/EtaXiLast;
  else
    gamma = 0.0;
  EtaXiLast = etaxi;
  
  Phi = Etap + gamma * LastPhi;
  LastPhi = Phi;
  
  Phipp = Phi - conjdot(c, Phi) * c;
  Phip = Phipp;
  Normalize (Phip);

  cerr << "E = " << E0 << endl;
}


void ConjGrad::Iterate()
{
  if (!IsSetup)
    Setup();
  // First, calculate steepest descent vector:
  //CalcPhiCG();
  CalcPhiSD();

  // Now, pick optimal theta for 
  // cnext = cos(theta)*c + sin(theta)*Phi;
  double dE_dtheta = 2.0*realconjdot(Phip, Hc);
  double theta1 = M_PI/300.0;
  cnext = cos(theta1)*c + sin(theta1)*Phip;
  H.Apply (cnext, Hcnext);
  double E1 = realconjdot (cnext, Hcnext);
  double A1 = (E0 - E1 + 0.5*sin(2.0*theta1)*dE_dtheta)/(1.0-cos(2.0*theta1));
  double B1 = 0.5*dE_dtheta;
  double thetaMin = 0.5*atan (B1/A1);
  //  cerr << "thetaMin = " << thetaMin << endl;

  c = cos(thetaMin)*c + sin(thetaMin)*Phip;
}
