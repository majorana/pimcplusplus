#include "NodalActionClass.h"
#include "../PathDataClass.h"
#include "../Common/MatrixOps/MatrixOps.h"


double 
FreeNodalActionClass::ActionImageSum (double L, double lambdaBeta, 
				      double disp)
{
  int numImages = 10;
  double sum = 0.0;
  double fourLambdaBetaInv = (lambdaBeta!=0.0) ?  1.0/(4.0*lambdaBeta) : 0.0;
  for (int image=-numImages; image<numImages; image++) {
    double x = disp + (double)image*L;
    sum += exp (-(x*x)*fourLambdaBetaInv);
  }
  return (-log(sum));
}

void FreeNodalActionClass::SetupFreeActions()
{
  const int nPoints = 1000;
  // Setup grids
  for (int i=0; i<NDIM; i++)
    ActionGrids[i].Init (-0.5*Path.GetBox()[i], 0.5*Path.GetBox()[i], nPoints);

  /// DEBUG
  FILE *fout = fopen ("FreeActions.dat", "w");

  Array<double,1> actionData(nPoints);
  // Now, setup up actions
  int nSplines = Path.TotalNumSlices/2 + (Path.TotalNumSlices%2)+1;
  ActionSplines.resize(nSplines);
  double lambdaTau = Path.tau * Path.Species(SpeciesNum).lambda;
  for (int spline=0; spline<nSplines; spline++) {
    double lambdaBeta = lambdaTau * (double)spline;
    for (int dim=0; dim<NDIM; dim++) {
      double L = Path.GetBox()[dim];
      for (int i=0; i<nPoints; i++) {
	double disp = ActionGrids[dim](i);
	actionData(i) = ActionImageSum (L, lambdaBeta, disp);
	fprintf (fout, "%1.12e ", actionData(i));
      }
      fprintf (fout, "\n");
      // Since the action is periodic, the slope should be zero
      // at the boundaries
      ActionSplines(spline)[dim].Init (&ActionGrids[dim], actionData,
				       0.0, 0.0);
    }
  }
  fclose (fout);
}



FreeNodalActionClass::FreeNodalActionClass (PathDataClass &pathData,
					    int speciesNum) :
  ActionBaseClass (pathData), 
  Path (pathData.Path), 
  SpeciesNum (speciesNum)
{
  int N = Path.Species(speciesNum).LastPtcl - 
      Path.Species(speciesNum).FirstPtcl+1;
  DetMatrix.resize(N,N);
  Cofactors.resize(N,N);
  GradVec.resize(N);
}


void 
FreeNodalActionClass::GradientDet (int slice, double &det, 
				   Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  int myStartSlice, myEndSlice;
  int myProc = PathData.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  double t = abs(refSlice-slice) * PathData.Action.tau;
  double beta = PathData.Path.TotalNumSlices * PathData.Action.tau;
  int sliceDiff = abs(slice-refSlice);
  sliceDiff = min (sliceDiff, Path.TotalNumSlices-sliceDiff);
  assert (sliceDiff <= Path.TotalNumSlices);
  assert (sliceDiff > 0);
  t = min (t, fabs(beta-t));
  assert (t <= 0.500000001*beta);
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  // HACK HACK HACK for now;  should work for serial mode.
  //   if (Path.GetRefSlice() < Path.NumTimeSlices())
  //     for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
  //       Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);
  //   Path.RefPath.AcceptCopy();

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      double action = 0.0;
      for (int dim=0; dim<NDIM; dim++)
	action += ActionSplines(sliceDiff)[dim](diff[dim]);
      //DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dist*dist);
      DetMatrix(refPtcl-first, ptcl-first) = exp(-action);
    }
  }


//   cerr << "slice = " << slice << endl;
//   cerr << "RefSlice = " << Path.GetRefSlice() << endl;
//   cerr << "DetMatrix = " << endl << DetMatrix << endl;


  // Compute determinant
  det = Determinant (DetMatrix);
  Cofactors = DetMatrix;
  GJInverse (Cofactors);
  Transpose (Cofactors);
  Cofactors = det * Cofactors;

//   double cofDet = 0.0;
//   for (int col=0; col<DetMatrix.cols(); col++)
//     cofDet += DetMatrix(col,0) * Cofactors(col,0);

  // Now compute gradient of determinant
  for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
    gradient(ptcl-first) = 0.0;
    dVec &r = Path(slice, ptcl);
    for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      dVec gradPhi;
      for (int dim=0; dim<NDIM; dim++)
	gradPhi[dim] = -ActionSplines(sliceDiff)[dim].Deriv(diff[dim]) 
	  * DetMatrix(refPtcl-first, ptcl-first);
      //dVec gradPhi = -2.0*C*diff*DetMatrix(refPtcl-first,ptcl-first);

      gradient(ptcl-first) = 
	gradient(ptcl-first)+ gradPhi*Cofactors(refPtcl-first, ptcl-first);
    }
  }
}



void 
FreeNodalActionClass::GradientDetFD (int slice, double &det, 
				     Array<dVec,1> &gradient)
{
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  // Fill up determinant matrix
  int myStartSlice, myEndSlice;
  int myProc = PathData.Communicator.MyProc();
  Path.SliceRange (myProc, myStartSlice, myEndSlice);
  int refSlice = Path.GetRefSlice()-myStartSlice;
  double t = abs(refSlice-slice) * PathData.Action.tau;
  double beta = PathData.Path.TotalNumSlices * PathData.Action.tau;
  t = min (t, fabs(beta-t));
  double lambda = species.lambda;
  double C = 1.0/(4.0*M_PI * lambda * t);

  // HACK HACK HACK for now;  should work for serial mode.
//   if (Path.GetRefSlice() < Path.NumTimeSlices())
//     for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++)
//       Path.RefPath(ptcl) = Path(Path.GetRefSlice(), ptcl);

  for (int refPtcl=species.FirstPtcl; refPtcl<=species.LastPtcl; refPtcl++) {
    const dVec &rRef = Path.RefPath(refPtcl);
    for (int ptcl=species.FirstPtcl; ptcl<=species.LastPtcl; ptcl++) {
      const dVec &r = Path(slice, ptcl);
      //      dVec diff = r-rRef;
      dVec diff;
      double dist;
      Path.RefDistDisp (slice, refPtcl, ptcl, dist, diff);
      DetMatrix(refPtcl-first, ptcl-first) = exp(-C*dist*dist);
    }
  }


  // Compute determinant
  det = Determinant (DetMatrix);

  dVec disp;
  double dist;

  double eps = 1.0e-5;
  for (int ptcl=first; ptcl <= last; ptcl++) {
    dVec delta = 0.0;
    dVec &r = Path(slice, ptcl);
    double dplus, dminus;
    for (int dim=0; dim<NDIM; dim++) {
      delta = 0.0;
      delta[dim] = 0.5*eps;
      for (int ref=first; ref <= last; ref++) {
	Path.RefDistDisp (slice, ref, ptcl, dist, disp);
	dVec diff = disp + delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dplus = Determinant (DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	Path.RefDistDisp (slice, ref, ptcl, dist, disp);
	dVec diff = disp + delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      dminus = Determinant(DetMatrix);
      for (int ref=first; ref <= last; ref++) {
	Path.RefDistDisp (slice, ref, ptcl, dist, disp);
	dVec diff = disp + delta;
        DetMatrix(ref-first, ptcl-first) = exp(-C*dot(diff,diff));
      }
      gradient(ptcl-first)[dim] = (dplus-dminus)/eps;
    }
  }
}


double FreeNodalActionClass::NodalDist (int slice)
{
  double grad2 = 0.0;
  SpeciesClass &species = Path.Species(SpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  double det;
  int N = last-first+1;

  GradientDet (slice, det, GradVec);
  //Array<dVec,1> gradFD(N);
  //GradientDetFD (slice, det, gradFD);
//   for (int i=0; i<N; i++)
//     for (int dim=0; dim<NDIM; dim++)
//       assert (fabs(gradFD(i)[dim] - GradVec(i)[dim]) < 1.0e-9);
//     fprintf (stderr, "(%1.6e %1.6e %1.6e) (%1.6e %1.6e %1.6e) \n", 
// 	     gradFD(i)[0], gradFD(i)[1], gradFD(i)[2], 
// 	     GradVec(i)[0], GradVec(i)[1], GradVec(i)[2]);
  for (int i=0; i<N; i++)
    grad2 += dot (GradVec(i), GradVec(i));
  return (det/sqrt(grad2));
}


void FreeNodalActionClass::Read (IOSectionClass &in)
{
  // Do nothing for now
  SetupFreeActions();
}



double FreeNodalActionClass::Action (int startSlice, int endSlice,
				     const Array<int,1> &changePtcls, 
				     int level)
{ 
  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Action.tau * (double)skip;

  int myStart, myEnd;
  Path.SliceRange(PathData.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;
  
  double dist1, dist2;
  if (startSlice != refSlice) {
    dist1 = NodalDist (startSlice);
    if (dist1 < 0.0)
      return 1.0e100;
  }
  else
    dist1 = sqrt(-1.0);
  
  int totalSlices = Path.TotalNumSlices;
  double uNode=0.0;
  for (int slice=startSlice; slice < endSlice; slice+=skip) {
    if ((slice+skip == refSlice) || (slice+skip == refSlice+totalSlices))
      dist2 = sqrt(-1.0);
    else {
      dist2 = NodalDist (slice+skip);
      if (dist2 < 0.0)
	return 1.0e100;
    }

    if (isnan (dist1) || (dist1==0.0))
      uNode -= log1p(-exp(-dist2*dist2/(lambda*levelTau)));
    else if (isnan(dist2) || (dist2==0.0))
      uNode -= log1p(-exp(-dist1*dist1/(lambda*levelTau)));
    else
      uNode -= log1p(-exp(-dist1*dist2/(lambda*levelTau)));
    dist1 = dist2;
  }
  return uNode;
}



double FreeNodalActionClass::d_dBeta (int slice1, int slice2, int level)
{ 
  SpeciesClass &species = Path.Species(SpeciesNum);
  double lambda = species.lambda;
  int skip = 1<<level;
  double levelTau = PathData.Action.tau * (double)skip;

  int myStart, myEnd;
  Path.SliceRange(PathData.Communicator.MyProc(), myStart, myEnd);
  int refSlice = Path.GetRefSlice() - myStart;
  
  double dist1, dist2;
  if (slice1 != refSlice) {
    dist1 = NodalDist (slice1);
    if (dist1 < 0.0)
      return 1.0e100;
  }
  else
    dist1 = sqrt(-1.0);
  
  int totalSlices = Path.TotalNumSlices;
  double uNode=0.0;
  for (int slice=slice1; slice < slice2; slice+=skip) {
    if ((slice+skip == refSlice) || (slice+skip == refSlice+totalSlices))
      dist2 = sqrt(-1.0);
    else {
      dist2 = NodalDist (slice+skip);
      if (dist2 < 0.0)
	return 1.0e100;
    }

    double prod;
    if (isnan (dist1) || (dist1==0.0))
      prod = dist2*dist2;
    else if (isnan(dist2) || (dist2==0.0))
      prod = dist1*dist1;
    else
      prod = dist1*dist2;

    uNode += prod/(lambda*levelTau*levelTau)/expm1(prod/(lambda*levelTau));
    dist1 = dist2;
  }
  return uNode/(double)Path.TotalNumSlices;
}
