#include "LangevinMove.h"

void
LangevinMoveClass::AccumForces()
{
  PathData.Actions.GetForces(Particles, FShortSum, FLongSum);
}


void
LangevinMoveClass::LDStep()
{
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

  assert(in.ReadVar("Mass", Mass));
  assert (in.ReadVar("TimeStep", TimeStep));
  assert (in.ReadVar("NumEquilSteps", NumEquilSteps));
  assert (in.ReadVar("NumAccumSteps", NumAccumSteps));
  assert (in.ReadVar("Species", speciesStr));

  LDSpecies = PathData.Path.SpeciesNum(speciesStr);
  assert (LDSpecies != -1);
  SpeciesClass &species = PathData.Path.Species(LDSpecies);
  int numPtcls = species.NumParticles;
  V.resize(numPtcls);
  R.resize(numPtcls);
  FShortSum.resize(numPtcls);
  FLongSum.resize(numPtcls);
  FShortTmp.resize(numPtcls);
  FLongTmp.resize(numPtcls);
  OldFShort.resize(numPtcls);
  OldFLong.resize(numPtcls);
  Particles.resize(numPtcls);
  WriteArray.resize(numPtcls,NDIM);

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
  Vsum = Vsum/(double)V.size();
  for (int i=0; i<V.size(); i++) {
    V(i) -= Vsum;
    Esum    += 0.5 * Mass * dot(V(i), V(i));
  }

  /// Now, normalize to appropriate temperature using equipartition
  /// theorem 
  double norm = sqrt(0.5*(double)(NDIM*V.size())*kBT/Esum);
  for (int i=0; i<V.size(); i++)
    V(i) = norm * V(i);

  /// Double-check that we did this write
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
