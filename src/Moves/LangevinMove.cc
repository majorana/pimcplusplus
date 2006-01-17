#include "LangevinMove.h"
#include <Common/MatrixOps/MatrixOps.h>


/// AccumForces adds the forces for my section of the path to the
/// deque of forces.
void
LangevinMoveClass::AccumForces()
{
  FShort = dVec(0.0);
  FLong  = dVec(0.0);
  PathData.Actions.GetForces(Particles, FShort, FLong);
  FShortSum += FShort;
  FLongSum  += FLong;
  /// Queue the forces for the computation of the variance.  Note that
  /// the copy is very important since the normal blitz copy constructor
  /// operator references the array data.
  FTmp += FShort + FLong;
  FDeque.push_back(FTmp.copy());
}

void 
Stats (Array<double,1> &x, double &mean, double &var, double &kappa)
{
  int N = x.size();
  double ninv = 1.0/(N);
  double xsum, x2sum;
  xsum = x2sum = 0.0;
  for (int i=0; i < N; i++) {
    xsum  += x(i);
    x2sum += x(i)*x(i);
  }
  mean = xsum * ninv;
  var = x2sum/(double)(N-1) - mean*mean;

  /// Compute autocorrelation time, kappa
  kappa = 1.0;
  bool done = false;
  int i = 1;
  while (!done) {
    double c_i = 0.0;
    for (int j=0; j<(N-i); j++)
      c_i += (x(j)-mean) * (x(j+i)-mean);
    c_i /= (var*(double)(N-i));
    i++;
    kappa += c_i;
    done = (i>=N) || (c_i<=0.0);
  }
}



void
LangevinMoveClass::CalcCovariance()
{
//   if (FDeque.size() < 20) {
//     A = 0.0;
//     for (int i=0; i<A.rows(); i++)
//       A(i,i) = 1.0e-5;
//     return;
//   }
  int N = R.size();
  double Ninv = 1.0/(double)N;

  /// First, sum the forces over the processors in my clone.
  for (int i=0; i<FDeque.size(); i++) {
    PathData.IntraComm.AllSum(FDeque[i], FTmp);
    FDeque[i] = FShort;
  }

  Fmean = 0.0;
  /// Compute mean force
  int n=0;
  for (int k=0; k<FDeque.size(); k++)
    for (int ptcl=0; ptcl<N; ptcl++) 
      for (int dim=0; dim<NDIM; dim++) {
	Fmean(n) += FDeque[k](ptcl)[dim];
	n++;
      }
  Fmean = (1.0/(double)FDeque.size())*Fmean;

  
  CoVar = 0.0;
  /// Now, calculate the covariance and 
  for (int k=0; k<FDeque.size(); k++) {
    int i=0;
    for (int ptcl1=0; ptcl1<N; ptcl1++) 
      for (int dim1=0; dim1<NDIM; dim1++) {
	int j=0;
	for (int ptcl2=0; ptcl2<N; ptcl2++)
	  for (int dim2=0; dim2<NDIM; dim2++) {
	    CoVar(i,j) += 
	      (FDeque[k](ptcl1)[dim1] * FDeque[k](ptcl2)[dim2]);
	    j++;
	  }
	i++;
      }
  }
  CoVar = (1.0/(double)(FDeque.size()-1))*CoVar;
  /// Subtract <F(i)><F(j)>
  for (int i=0; i<CoVar.rows(); i++)
    for (int j=0; j<CoVar.cols(); j++)
      CoVar(i,j) -= Fmean(i)*Fmean(j);

  /// Compute the average autocorrelation time for the diagaonal
  /// elements only.
  Array<double,1> x(FDeque.size());
  double tmpMean, tmpVar, tmpKappa;
  kappa = 0.0;
  for (int ptcl=0; ptcl<N; ptcl++)
    for (int dim=0; dim<NDIM; dim++) {
      for (int k=0; k<x.size(); k++) 
	x(k) = FDeque[k](ptcl)(dim);
      Stats (x, tmpMean, tmpVar, tmpKappa);
      kappa += tmpKappa;
    }
  kappa /= (double)(N*NDIM);

  if (PathData.IntraComm.MyProc()==0)
    cerr << "CloneNum = " << PathData.GetCloneNum()
	 << " and kappa = " << kappa << endl;

  /// Multiply covariance by kappa
  CoVar *= kappa;

  /// Now we have the covariance matrix for the clone.  We must sum
  /// over all the clones and divide by the number
  if (PathData.IntraComm.MyProc() == 0) {
    int numClones = PathData.GetNumClones();

    PathData.InterComm.AllSum(Fmean, FallSum);
    Fmean = 1.0/(double)numClones *FallSum;

    PathData.InterComm.AllSum(CoVar, A);
    double beta = PathData.Path.tau * PathData.Path.TotalNumSlices;
    A *= 0.5*TimeStep*beta/(Mass*(double)(numClones*numClones));
    /// Now calculate eigenvalues and eigenvectors
    SymmEigenPairs (A, A.extent(0), Lambda, L);
  }
  /// Now, broadcast eigenvectors and values to all processors
  PathData.IntraComm.Broadcast(0, L);
  PathData.IntraComm.Broadcast(0, Lambda);
  PathData.IntraComm.Broadcast(0, Fmean);

  /// Transpose eigenvectors
  for (int i=0; i<L.rows(); i++)
    for (int j=0; j<L.cols(); j++)
      Ltrans(i,j) = L(j,i);

  /// Now empty out the force deque for next time
  FDeque.clear();

}




/// Friction sets
void
LangevinMoveClass::CalcFriction()
{
  if (NumFs < 20) {
    Lambda = 1.0e-5;
    L      = 0.0;
    Ltrans = 0.0;
    for (int i=0; i<L.extent(0); i++) {
      L(i,i) = 1.0;
      Ltrans(i,i) = 1.0;
    }
  }
  else {
    double ninv = 1.0/NumFs;
    double beta = PathData.Path.tau * PathData.Path.TotalNumSlices;
    double norm = 0.5*beta*Mass;
    for (int i=0; i<A.extent(0); i++)
      for (int j=0; j<A.extent(1); j++)
	A(i,j) = (ninv*(FF(i,j) - ninv*Fmean(i)*Fmean(j)));
    SymmEigenPairs(A, A.extent(0), Lambda, L);
    for (int i=0; i<L.extent(0); i++)
      for (int j=0; j<L.extent(1); i++)
	Ltrans(i,j) = L(j,i);
  }
}


void
LangevinMoveClass::LDStep()
{
  CalcCovariance();
  /// Write out positions and velocities
  TimeVar.Write(Time);
  for (int i=0; i<R.size(); i++) {
    dVec r = R(i);
    PathData.Path.PutInBox(r);
    for (int j=0; j<NDIM; j++) 
      WriteArray(i,j) = r[j];
  }
  Rvar.Write(WriteArray);

  for (int i=0; i<V.size(); i++) 
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = V(i)[j];
  Vvar.Write(WriteArray);

  for (int i=0; i<OldFShort.size(); i++) 
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = OldFShort(i)[j];
  FShortVar.Write(WriteArray);

  for (int i=0; i<OldFShort.size(); i++) 
    for (int j=0; j<NDIM; j++)
      WriteArray(i,j) = OldFLong(i)[j];
  FLongVar.Write(WriteArray);
  FLongVar.Flush();

  // OldF holds the force computed at x(t).
  for (int i=0; i<R.size(); i++)
    R(i) += TimeStep * V(i) + 
      0.5*TimeStep*TimeStep*MassInv*(OldFShort(i)+OldFLong(i));

  // Now, compute F(t+dt)
  // Sum forces over all the processors (all clones included)
  PathData.WorldComm.AllSum(FShortSum, FShortTmp);
  PathData.WorldComm.AllSum(FLongSum,  FLongTmp);
  int numClones = PathData.GetNumClones();
  double norm = 1.0/(double)(numClones*NumAccumSteps);
  for (int i=0; i<FShortSum.size(); i++) {
    FShortSum(i) = norm*FShortTmp(i);
    FLongSum(i) = norm*FLongTmp(i);
  }
  // Now FShortSum and FLongSum hold the forces at x(t+dt)
  
  // Compute V(t+dt)
  for (int i=0; i<V.size(); i++)
    V(i) += 0.5*MassInv*TimeStep*
      (FShortSum(i) + FLongSum(i) + OldFShort(i) + OldFLong(i));

  // Copy new force into old
  OldFShort = FShortSum; 
  OldFLong = FLongSum;

  // And reset FShortSum and FLongSum
  dVec zero(0.0);
  for (int i=0; i<FShortSum.size(); i++) {
    FShortSum(i) = zero;
    FLongSum(i) = zero;
  }

  /// Accumulate statistics for covariance matrices.
  for (int i=0; i<FF.extent(0); i++) {
    int m1 = i/NDIM;
    int m2 = i%NDIM;
    Fmean(i) += OldFShort(m1)[m2] + OldFLong(m1)[m2];
    for (int j=0; j<FF.extent(1); j++) {
      int n1 = j/NDIM;
      int n2 = j%NDIM;
      FF(i,j) += ( (OldFShort(m1)[m2] + OldFLong(m1)[m2]) * 
		   (OldFShort(n1)[n2] + OldFLong(n1)[n2]) );
    }
  }
  NumFs++;


  // Put x(t+2dt) into the Path so we can start accumulating forces
  // for the next step
  int first = PathData.Path.Species(LDSpecies).FirstPtcl;
  SetMode (NEWMODE);
  for (int slice=0; slice<PathData.Path.NumTimeSlices(); slice++)
    for (int i=0; i<R.size(); i++) 
      PathData.Path(slice,i+first) = R(i) + TimeStep*V(i) + 
	0.5*TimeStep*TimeStep* (OldFShort(i) + OldFLong(i));

  /// Warp electron paths to follow ions
  PathData.Path.WarpPaths(LDSpecies);
  
  SetMode (OLDMODE);
  for (int slice=0; slice<PathData.Path.NumTimeSlices(); slice++)
    for (int i=0; i<R.size(); i++) 
      PathData.Path(slice,i+first) = R(i) + TimeStep*V(i) + 
	0.5*TimeStep*TimeStep*(OldFShort(i)+OldFLong(i));

  /// Increment the time
  Time += TimeStep;
  
  /// Update the nodal actions
  PathData.Actions.UpdateNodalActions();
}

void 
LangevinMoveClass::MakeMove()
{
  if (MCSteps >= (NumEquilSteps+NumAccumSteps)) {
    MCSteps = 0;
    LDStep();
  }
  else if (MCSteps >= NumEquilSteps) {
    AccumForces();
    MCSteps++;
  }
  else
    MCSteps++;
}

void 
LangevinMoveClass::Read(IOSectionClass &in)
{
  string speciesStr;

  assert (in.ReadVar("Mass",          Mass));
  assert (in.ReadVar("TimeStep",      TimeStep));
  assert (in.ReadVar("NumEquilSteps", NumEquilSteps));
  assert (in.ReadVar("NumAccumSteps", NumAccumSteps));
  assert (in.ReadVar("Species",       speciesStr));

  LDSpecies = PathData.Path.SpeciesNum(speciesStr);
  assert (LDSpecies != -1);
  SpeciesClass &species = PathData.Path.Species(LDSpecies);
  int numPtcls = species.NumParticles;
  V.resize         (numPtcls);
  R.resize         (numPtcls);
  FShort.resize    (numPtcls);
  FLong.resize     (numPtcls);
  FTmp.resize      (numPtcls);
  FShortSum.resize (numPtcls);
  FLongSum.resize  (numPtcls);
  FShortTmp.resize (numPtcls);
  FLongTmp.resize  (numPtcls);
  OldFShort.resize (numPtcls);
  OldFLong.resize  (numPtcls);
  Particles.resize (numPtcls);
  WriteArray.resize(numPtcls,NDIM);
  Fmean.resize  (NDIM*numPtcls);
  FallSum.resize(NDIM*numPtcls);
  FF.resize     (NDIM*numPtcls, NDIM*numPtcls);
  A.resize      (NDIM*numPtcls, NDIM*numPtcls);
  CoVar.resize  (NDIM*numPtcls, NDIM*numPtcls);
  L.resize      (NDIM*numPtcls, NDIM*numPtcls);
  Ltrans.resize (NDIM*numPtcls, NDIM*numPtcls);
  Lambda.resize (NDIM*numPtcls);
  NumFs = 0;

  dVec zero(0.0);
  for (int i=0; i<numPtcls; i++) {
    Particles(i) = i + species.FirstPtcl;
    R(i) = PathData.Path(0,Particles(i));
    FShortSum(i) = zero;
    FLongSum(i)  = zero;
    OldFShort(i) = zero;
    OldFLong(i)  = zero;
  }
  MassInv = 1.0/Mass;


  /// Initialize the velocities
  InitVelocities();

}


void LangevinMoveClass::InitVelocities()
{
  double beta = PathData.Path.tau * PathData.Path.TotalNumSlices;
  double kBT = 1.0/beta;
  /// First, initialize randomly.  We must use global random numbers,
  /// since all clones must have the same velocities.
  dVec Vsum;
  for (int i=0; i<NDIM; i++)
    Vsum[i] = 0.0;
  for (int i=0; i<V.size(); i++) {
    PathData.Path.Random.WorldGaussianVec(1.0, V(i));
    Vsum += V(i);
  }
  // Now, remove center-of-mass velocity
  double Esum = 0.0;
  Vsum = (1.0/(double)V.size())*Vsum;
  for (int i=0; i<V.size(); i++) {
    V(i) -= Vsum;
    Esum    += 0.5 * Mass * dot(V(i), V(i));
  }

  /// Now, normalize to appropriate temperature using equipartition
  /// theorem 
  double norm = sqrt(0.5*(double)(NDIM*V.size())*kBT/Esum);
  for (int i=0; i<V.size(); i++)
    V(i) = norm * V(i);

  /// Double-check that we did this right
  for (int i=0; i<NDIM; i++)
    Vsum[i] = 0.0;
  Esum = 0.0;
  for (int i=0; i<V.size(); i++) {
    Esum += 0.5 * Mass * dot(V(i), V(i));
    Vsum += V(i);
  }
  cerr << "Esum = " << Esum << " 3/2N kB T = " <<
    (0.5*(double)(NDIM*V.size())*kBT) << endl;
  assert (fabs((Esum/(0.5*(double)(NDIM*V.size())*kBT))-1.0) < 1.0e-12);
  assert ((0.5*Mass*dot(Vsum,Vsum)) < 1.0e-10*kBT);
}
