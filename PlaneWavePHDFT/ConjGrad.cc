#include "ConjGrad.h"


void ConjGrad::Setup()
{
  c.resize (H.GVecs.size());
  cnext.resize (H.GVecs.size());
  Hcnext.resize (H.GVecs.size());
  Hc.resize (H.GVecs.size());
  SDvec.resize (H.GVecs.size());
  c = 1.0;
  Normalize (c);
  IsSetup = true;
}

void ConjGrad::CalcSD()
{
  H.Apply(c, Hc);
  E0 = realconjdot (c, Hc);
  cerr << "E = " << E0 << endl;
  SDvec = E0*c - Hc;
  Normalize (SDvec);
}

void ConjGrad::Iterate()
{
  if (!IsSetup)
    Setup();
  // First, calculate steepest descent vector:
  CalcSD();

  // Now, pick optimal theta for cnext = cos(theta)*c +
  // sin(theta)*SDvec;
  double dE_dtheta = 2.0*realconjdot(SDvec, Hc);
  double theta1 = M_PI/300.0;
  cnext = cos(theta1)*c + sin(theta1)*SDvec;
  H.Apply (cnext, Hcnext);
  double E1 = realconjdot (cnext, Hcnext);
  

  double A1 = (E0 - E1 + 0.5*dE_dtheta)/(1.0-cos(2.0*theta1));
  double B1 = 0.5*dE_dtheta;
  double thetaMin = 50.0*0.5*atan (B1/A1);
  //  cerr << "thetaMin = " << thetaMin << endl;

  c = cos (thetaMin)*c + sin(thetaMin)*SDvec;
  Normalize (c);
}
