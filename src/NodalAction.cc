#include "NodalAction.h"
#include "Common/MatrixOps/MatrixOps.h"

void 
FPNodalActionClass::GradientDet (int slice, double &det, 
				 Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  double t = slice * PathData.Action.tau;
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    const dVec &rRef = Path(0, refPtcl);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      dVec diff = r-rRef;
      DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dot(diff,diff));
    }
  }

  // Compute determinant
  det = Determinant (DetMatrix);
  Cofactors = DetMatrix;
  GJInverse (Cofactors);
  Transpose (Cofactors);
  Cofactors = det * Cofactors;

  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl) = 0.0, 0.0, 0.0;
    dVec &r = Path(slice, ptcl);
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec &rRef = Path(0, refPtcl);
      dVec diff =  r-rRef;
      dVec gradPhi = -2.0*C*diff*DetMatrix(refPtcl-first,ptcl-first);
      gradient(ptcl-first) += Cofactors(refPtcl-first, ptcl-first)*gradPhi;
    }
  }
}



void 
FPNodalActionClass::GradientDetFD (int slice, double &det, 
				   Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  double t = slice * PathData.Action.tau;
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    const dVec &rRef = Path(0, refPtcl);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      dVec diff = r-rRef;
      DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dot(diff,diff));
    }
  }

  // Compute determinant
  det = Determinant (DetMatrix);

  
  double eps = 1.0e-4;
  for (int ptcl=first; ptcl <= last; ptcl++) {
    dVec delta = 0.0;
    dVec &r = Path(slice, ptcl);
    double dplus, dminus;
    for (int dim=0; dim<NDIM; dim++) {
      delta[dim] = 0.5*eps;
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path(slice, ref);
	dVec diff = r - rRef + delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dplus = Determinant (DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path(slice, ref);
	dVec diff = r - rRef - delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dminus = Determinant(DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path(slice, ref);
	dVec diff = r - rRef;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      gradient(ptcl-first)[dim] = (dplus-dminus)/eps;
    }
  }
}


double FPNodalActionClass::NodalDist (int slice)
{
  double dist2;
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  double det;
  int N = last-first+1;
  Array<dVec,1> grad(N);
  GradientDet (slice, det, GradVec);
  for (int i=0; i<N; i++)
    dist2 += dot (GradVec(i), GradVec(i));
  return (det/sqrt(dist2));
}


double FPNodalActionClass::Action (int startSlice, int endSlice, int level,
				   Array<int,1> &changePtcls)
{ 
  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Action.tau * (double)skip;

  double dist1, dist2;
  dist1 = NodalDist (startSlice);
  if (dist1 < 0.0)
    return 1.0e100;
  double uNode=0.0;
  for (int slice=startSlice; slice < endSlice; slice++) {
    dist2 = NodalDist (slice+skip);
    if (dist2 < 0.0)
      return 1.0e100;
    uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
    dist1 = dist2;
  }
  
  return uNode;
}

