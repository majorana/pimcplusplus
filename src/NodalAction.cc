#include "NodalAction.h"
#include "Common/MatrixOps/MatrixOps.h"

void 
FPNodalActionClass::GradientDet (int slice, double &det, 
				 Array<dVec,1> &gradient)
{
  if ((slice == Path.GetRefSlice()) || (slice == (Path.GetRefSlice()+Path.TotalNumSlices))) {
    det = 1.0;
    gradient = 1.0e-10;
    return;
  }

  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  double t = abs(slice-PathData.Path.GetRefSlice()) * PathData.Action.tau;
  double beta = PathData.Path.TotalNumSlices * PathData.Action.tau;
  t = min (t, fabs(beta-t));
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  // HACK HACK HACK for now;  should work for serial mode.
  if (Path.GetRefSlice() < Path.NumTimeSlices())
    for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
      Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    const dVec &rRef = Path.RefPath(refPtcl);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      dVec diff = r-rRef;
      DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dot(diff,diff));
    }
  }

  cerr << "slice = " << slice << endl;
  cerr << "RefSlice = " << Path.GetRefSlice() << endl;
  cerr << "DetMatrix = " << endl << DetMatrix << endl;


  // Compute determinant
  det = Determinant (DetMatrix);
  Cofactors = DetMatrix;
  GJInverse (Cofactors);
  Transpose (Cofactors);
  Cofactors = det * Cofactors;

  double cofDet = 0.0;
  for (int col=0; col<DetMatrix.cols(); col++)
    cofDet += DetMatrix(col,0) * Cofactors(col,0);
  cerr << "Determinant      = " << det << endl;
  cerr << "Det by cofactors = " << cofDet << endl;

  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl-first) = 0.0;
    cerr << gradient(ptcl-first) << endl;
    dVec &r = Path(slice, ptcl);
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec &rRef = Path.RefPath(refPtcl);
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
  if ((slice == Path.GetRefSlice()) || (slice == (Path.GetRefSlice()+Path.TotalNumSlices))) {
    det = 1.0;
    gradient = 1.0e-10;
    return;
  }
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  double t = slice * PathData.Action.tau;
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  // HACK HACK HACK for now;  should work for serial mode.
  if (Path.GetRefSlice() < Path.NumTimeSlices())
    for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
      Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    const dVec &rRef = Path.RefPath(refPtcl);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      dVec diff = r-rRef;
      DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dot(diff,diff));
    }
  }


  // Compute determinant
  det = Determinant (DetMatrix);

  
  double eps = 1.0e-5;
  for (int ptcl=first; ptcl <= last; ptcl++) {
    dVec delta = 0.0;
    dVec &r = Path(slice, ptcl);
    double dplus, dminus;
    for (int dim=0; dim<NDIM; dim++) {
      delta = 0.0;
      delta[dim] = 0.5*eps;
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path.RefPath(ref);
	dVec diff = r - rRef + delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dplus = Determinant (DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path.RefPath(ref);
	dVec diff = r - rRef - delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dminus = Determinant(DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	dVec &rRef = Path.RefPath(ref);
	dVec diff = r - rRef;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      gradient(ptcl-first)[dim] = (dplus-dminus)/eps;
    }
  }
}


double FPNodalActionClass::NodalDist (int slice)
{
  double grad2;
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  double det;
  int N = last-first+1;

  GradientDet (slice, det, GradVec);
  Array<dVec,1> gradFD(N);
  GradientDetFD (slice, det, gradFD);
  for (int i=0; i<N; i++)
    fprintf (stderr, "(%1.6e %1.6e %1.6e) (%1.6e %1.6e %1.6e) \n", 
	     gradFD(i)[0], gradFD(i)[1], gradFD(i)[2], 
	     GradVec(i)[0], GradVec(i)[1], GradVec(i)[2]);
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec(i), GradVec(i));
  return (det/sqrt(grad2));
}


double FPNodalActionClass::Action (int startSlice, int endSlice,
				   const Array<int,1> &changePtcls, int level)
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
  for (int slice=startSlice; slice < endSlice; slice+=skip) {
    dist2 = NodalDist (slice+skip);
    if (dist2 < 0.0)
      return 1.0e100;
    uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
    dist1 = dist2;
  }
  
  return uNode;
}

