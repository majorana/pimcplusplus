/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

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
  FTmp = FShort + FLong;
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
  var  = ninv * x2sum - mean*mean;
  var *= (double)N/(double)(N-1);

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
    if (c_i > 0)
      kappa += 2.0*c_i;
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
  int numFs = FDeque.size();
  double Ninv = 1.0/(double)N;
  
  /// First, sum the forces over the processors in my clone.
  for (int i=0; i<numFs; i++) {
    FTmp = dVec(0.0);
    PathData.IntraComm.AllSum(FDeque[i], FTmp);
    FDeque[i] = FTmp;
  }
  
  Fmean = 0.0;
  // Compute mean force
  for (int k=0; k<numFs; k++) {
    int n=0;
    for (int ptcl=0; ptcl<N; ptcl++) 
      for (int dim=0; dim<NDIM; dim++) {
	Fmean(n) += FDeque[k](ptcl)[dim];
	n++;
      }
  }
  Fmean = (1.0/(double)numFs)*Fmean;

  
  CoVar = 0.0;
  /// Now, calculate the covariance and 
  for (int k=0; k<numFs; k++) {
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
  CoVar = (1.0/(double)numFs)*CoVar;
  /// Subtract <F(i)><F(j)>
  for (int i=0; i<CoVar.rows(); i++)
    for (int j=0; j<CoVar.cols(); j++)
      CoVar(i,j) -= Fmean(i)*Fmean(j);
  /// Necessary to correct for bias that comes from the fact that we
  /// calculated the mean from the data.
  CoVar = ((double)numFs/(double)(numFs-1)) * CoVar;

  /// Compute the average autocorrelation time for the diagaonal
  /// elements only.
    Array<double,1> x(numFs);
  double tmpMean, tmpVar, tmpKappa;
  kappa = 0.0;
  for (int ptcl=0; ptcl<N; ptcl++)
    for (int dim=0; dim<NDIM; dim++) {
      for (int k=0; k<x.size(); k++) 
	x(k) = FDeque[k](ptcl)[dim];
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
    /// Make sure we don't have negative eigenvalues
    // cerr << "Lambdas = [" << Lambda << "]\n";
    //    cerr << "Fmean = " << Fmean << endl;
    for (int i=0; i<Lambda.size(); i++) {
      if (Lambda(i) <= 0.0) {
	cerr << "Negative eigenvalue detected. Lambda = "
	     << Lambda(i) << endl;
	Lambda(i) = 1.0e-10;
      }
    }
    /////////////////////////
    // HACK HACK HACK HACK //
    /////////////////////////
//     Lambda = 1.0e-10;
//     L = 0.0;
//     for (int i=0; i<L.rows(); i++)
//       L(i,i) = 1.0;
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




void
LangevinMoveClass::CalcCovariance2()
{
  int N = R.size();
  int numFs = FDeque.size();
  double Ninv = 1.0/(double)N;
  
  /// First, sum the forces over the processors in my clone.
  for (int i=0; i<numFs; i++) {
    FTmp = dVec(0.0);
    PathData.IntraComm.AllSum(FDeque[i], FTmp);
    FDeque[i] = FTmp;
  }
  
  Fmean = 0.0;
  // Compute mean force
  for (int k=0; k<numFs; k++) {
    int n=0;
    for (int ptcl=0; ptcl<N; ptcl++) 
      for (int dim=0; dim<NDIM; dim++) {
	Fmean(n) += FDeque[k](ptcl)[dim];
	n++;
      }
  }
  Fmean = (1.0/(double)numFs)*Fmean;


  if (PathData.IntraComm.MyProc() == 0) {
    int numC = PathData.GetNumClones();
    Array<double,2> meanMat(numC, Fmean.size());
    for (int k=0; k<Fmean.size(); k++)
      meanMat(PathData.GetCloneNum(), k) = Fmean(k);
    
    PathData.InterComm.AllGatherRows(meanMat);
    CoVar = 0.0;
    Fmean = 0.0;
    for (int k=0; k<meanMat.rows(); k++)
      Fmean += meanMat(k,Range::all());
    Fmean *= (1.0/(double)numC);

    for (int i=0; i<meanMat.cols(); i++)
      for (int j=0; j<meanMat.cols(); j++)
	for (int k=0; k<meanMat.rows(); k++)
	  CoVar(i,j) += meanMat(k,i)*meanMat(k,j);
    
    CoVar = (1.0/(double)numC)*CoVar;
    for (int i=0; i<CoVar.rows(); i++)
      for (int j=0; j<CoVar.cols(); j++)
	CoVar(i,j) -= Fmean(i)*Fmean(j);
    // Want Covar to be the square of the error.  Therefore, divide by
    // N-1
    CoVar = (1.0/(double)(numC-1)) * CoVar;
    // Compensatory term for "extra" noise
    for (int i=0; i<CoVar.rows(); i++)
      CoVar(i,i) += ExtraNoiseSigma*ExtraNoiseSigma;
    double beta = PathData.Path.tau * PathData.Path.TotalNumSlices;
    A = (0.5 * TimeStep*beta/Mass)*CoVar;
    // Now calculate eigenvalues and eigenvectors
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
  
  for (int i=0; i<Fmean.size(); i++)
    Fmean(i) += PathData.Path.Random.WorldGaussian(ExtraNoiseSigma);

  /// Now empty out the force deque for next time
  FDeque.clear();
}




/// Friction sets
void
LangevinMoveClass::CalcFriction()
{

}


void
LangevinMoveClass::VerletStep()
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
  Vec2Array (V, WriteArray);        Vvar.Write(WriteArray);
  Vec2Array(OldFShort, WriteArray); FShortVar.Write(WriteArray);
  Vec2Array (OldFLong, WriteArray); FLongVar.Write(WriteArray);
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
  /// And the Rho_ks
  PathData.Path.UpdateRho_ks();
}


void
MakePeriodic(Array<double,3> &A)
{
  int Nx = A.extent(0);
  int Ny = A.extent(1);
  int Nz = A.extent(2);
  for (int iy=0; iy<Ny-1; iy++)
    for (int iz=0; iz<Nz-1; iz++)
      A(Nx-1,iy,iz) = A(0,iy,iz);
  // Y face
  for (int ix=0; ix<Nx-1; ix++)
    for (int iz=0; iz<Nz-1; iz++)
      A(ix,Ny-1,iz) = A(ix,0,iz);
  // Z face
  for (int ix=0; ix<Nx-1; ix++)
    for (int iy=0; iy<Ny-1; iy++)
      A(ix,iy,Nz-1) = A(ix,iy,0);
  // XY edge
  for (int iz=0; iz<Nz-1; iz++)
    A(Nx-1,Ny-1,iz) = A(0,0,iz);
  // XZ edge
  for (int iy=0; iy<Ny-1; iy++)
    A(Nx-1,iy,Nz-1) = A(0,iy,0);
  // YZ edge
  for (int ix=0; ix<Nx-1; ix++)
    A(ix,Ny-1,Nz-1) = A(ix,0,0);
  /// Corner
  A(Nx-1,Ny-1,Nz-1) = A(0,0,0);
}

void
LangevinMoveClass::LangevinStep()
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
  Vec2Array (V, WriteArray);    Vvar.Write(WriteArray);
  Vec2Array (Vold, WriteArray); VOldVar.Write(WriteArray);
  CoVarVar.Write(CoVar);
  if (PathData.Actions.NodalActions(0) != NULL)
    if (PathData.Actions.NodalActions(0)->Type() == GROUND_STATE_FP) {
      FixedPhaseActionClass &FP = 
	*((FixedPhaseActionClass*) PathData.Actions.NodalActions(0));
      Array<double,1> bandEnergies;
  
      FP.GetBandEnergies(bandEnergies);
      BandEnergiesVar.Write(bandEnergies);
      if (DumpRho && (PathData.GetCloneNum()==0)) {
// 	FP.CalcDensity(Rho);
// 	RhoVar.Write(Rho);
	const Array<double,3> &rho = FP.GetDensity();
	int nx = rho.extent(0);
	int ny = rho.extent(1);
	int nz = rho.extent(2);
	Rho.resize(nx+1,ny+1, nz+1);
	Rho(Range(0,nx-1), Range(0,ny-1), Range(0,nz-1)) = rho;
	MakePeriodic(Rho);
	RhoVar.Write(Rho);
      }      
      if (DumpBandRho && (PathData.GetCloneNum()==0)) {
	FP.CalcBandDensity(BandRho);
	BandRhoVar.Write(BandRho);
      }
    }
  // Calculate the mean force, covariance and eigenvalue
  // decomposition. 
	
  if (PathData.GetNumClones() > 2)
    CalcCovariance2();
  else
    CalcCovariance();
  LambdaVar.Write(Lambda);
  Lambda = FrictionFactor * Lambda;
  for (int i=0; i<WriteArray.extent(0); i++) {
    WriteArray(i,0) = Fmean(3*i+0);
    WriteArray(i,1) = Fmean(3*i+1);
    WriteArray(i,2) = Fmean(3*i+2);
  }
  FVar.Write(WriteArray);
  FVar.Flush();


  for (int i=0; i<Lambda.size(); i++) {
    ExpLambda(i)   =   exp(-TimeStep*Lambda(i));
    Expm1Lambda(i) = expm1(-TimeStep*Lambda(i));
  }
  
  Vec2Array(R,  RArray);
  Vec2Array(Rp, RpArray);
  /// Store Fmean for later
  Array2Vec(Fmean, FTmp);

  // Compute next R;  from eq. 5.22 of attaccalite thesis
  /// Change into eigenbasis of covariace
  MatVecProd (Ltrans, RArray,  TempArray);  RArray  = TempArray;
  MatVecProd (Ltrans, RpArray, TempArray);  RpArray = TempArray;
  MatVecProd (Ltrans, Fmean,   TempArray);  Fmean   = TempArray;
  
  for (int i=0; i<RArray.size(); i++) {
    RArray (i) *= (1.0+ExpLambda(i));
    RpArray(i) *= ExpLambda(i);
    Fmean  (i) *= -MassInv*TimeStep/Lambda(i)*Expm1Lambda(i);
  }

  /// Change back into regular position space
  MatVecProd (L, RArray,  TempArray);  RArray  = TempArray;
  MatVecProd (L, RpArray, TempArray);  RpArray = TempArray;
  MatVecProd (L, Fmean,   TempArray);  Fmean   = TempArray;

  /// Update to new positions
  Rp = R;
  for (int i=0; i<R.size(); i++) {
    R(i)[0] = RArray(3*i+0) - RpArray(3*i+0) + Fmean(3*i+0);
    R(i)[1] = RArray(3*i+1) - RpArray(3*i+1) + Fmean(3*i+1);
    R(i)[2] = RArray(3*i+2) - RpArray(3*i+2) + Fmean(3*i+2);
  }

  /// Calculate velocity with 5.23 of Attacalite
  // First, do the R-Rp part
  Vec2Array(R,  RArray);
  Vec2Array(Rp,  RpArray);
  RArray -= RpArray;
  MatVecProd (Ltrans, RArray,  TempArray);  RArray  = TempArray;
  for (int i=0; i<RArray.size(); i++) 
    RArray (i) *= (-Lambda(i)*ExpLambda(i)/Expm1Lambda(i));
  MatVecProd (L, RArray,  TempArray);  RArray  = TempArray;
  /// Now do the force part
  Vec2Array(FTmp, Fmean);
  MatVecProd (Ltrans, Fmean, TempArray);    Fmean = TempArray;
  for (int i=0; i<Fmean.size(); i++) 
    Fmean(i) *= TimeStep*MassInv*
      (1.0 - (Expm1Lambda(i) + Lambda(i)*TimeStep)/
       (Lambda(i)*TimeStep*(-Expm1Lambda(i))));
  MatVecProd(L, Fmean, TempArray);      Fmean = TempArray;
  for (int i=0; i<V.size(); i++) {
    V(i)[0] = RArray(3*i+0) + Fmean(3*i+0);
    V(i)[1] = RArray(3*i+1) + Fmean(3*i+1);
    V(i)[2] = RArray(3*i+2) + Fmean(3*i+2);
  }
  /// This formula is inaccurate.  Use 5.23 of Attaccalite thesis
  Vold = (1.0/TimeStep)*(R-Rp);

  // Put x(t+2dt) into the Path so we can start accumulating forces
  // for the next step
  int first = PathData.Path.Species(LDSpecies).FirstPtcl;
  SetMode (NEWMODE);
  for (int slice=0; slice<PathData.Path.NumTimeSlices(); slice++)
    for (int i=0; i<R.size(); i++) 
      PathData.Path(slice,i+first) = R(i);
  // + TimeStep*V(i) + 0.5*TimeStep*TimeStep* (OldFShort(i) + OldFLong(i));

  /// Warp electron paths to follow ions
  PathData.Path.WarpPaths(LDSpecies);
  
  SetMode (OLDMODE);
  for (int slice=0; slice<PathData.Path.NumTimeSlices(); slice++)
    for (int i=0; i<R.size(); i++) 
      PathData.Path(slice,i+first) = R(i);
  // + TimeStep*V(i) + 0.5*TimeStep*TimeStep*(OldFShort(i)+OldFLong(i));

  /// Increment the time
  Time += TimeStep;
  
  /// Update the nodal actions
  PathData.Actions.UpdateNodalActions();

  /// And the Rho_ks
  PathData.Path.UpdateRho_ks();
}




void 
LangevinMoveClass::MakeMove()
{
  if (MCSteps >= (NumEquilSteps+NumAccumSteps)) {
    MCSteps = 0;
    if (Integrator == VERLET)
      VerletStep();
    else
      LangevinStep();
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
  if (in.ReadVar("ExtraNoiseSigma", ExtraNoiseSigma))
    perr << "Extra Langevin noise = " << ExtraNoiseSigma << endl;
  assert (in.ReadVar("TimeStep",      TimeStep));
  assert (in.ReadVar("NumEquilSteps", NumEquilSteps));
  assert (in.ReadVar("NumAccumSteps", NumAccumSteps));
  assert (in.ReadVar("Species",       speciesStr));
  if (in.ReadVar("FrictionFactor", FrictionFactor))
    perr << "FrictionFactor = " << FrictionFactor << endl;
  string integrator;
  bool found = in.ReadVar ("Integrator", integrator);
  if (found) {
    if (integrator == "VERLET")
      Integrator = VERLET;
    else if (integrator == "LANGEVIN")
      Integrator = LANGEVIN;
    else {
      cerr << "Unrecognized MD integrator \"" << integrator 
	   << "\".  Exitting.\n";
      abort();
    }
  }

  LDSpecies = PathData.Path.SpeciesNum(speciesStr);
  assert (LDSpecies != -1);
  SpeciesClass &species = PathData.Path.Species(LDSpecies);
  int numPtcls = species.NumParticles;
  V.resize           (numPtcls);
  Vold.resize        (numPtcls);
  R.resize           (numPtcls);
  Rp.resize          (numPtcls);
  FShort.resize      (numPtcls);
  FLong.resize       (numPtcls);
  FTmp.resize        (numPtcls);
  FShortSum.resize   (numPtcls);
  FLongSum.resize    (numPtcls);
  FShortTmp.resize   (numPtcls);
  FLongTmp.resize    (numPtcls);
  OldFShort.resize   (numPtcls);
  OldFLong.resize    (numPtcls);
  Particles.resize   (numPtcls);
  WriteArray.resize  (numPtcls,NDIM);
  Fmean.resize       (NDIM*numPtcls);
  FallSum.resize     (NDIM*numPtcls);
  A.resize           (NDIM*numPtcls, NDIM*numPtcls);
  CoVar.resize       (NDIM*numPtcls, NDIM*numPtcls);
  L.resize           (NDIM*numPtcls, NDIM*numPtcls);
  Ltrans.resize      (NDIM*numPtcls, NDIM*numPtcls);
  Lambda.resize      (NDIM*numPtcls);
  ExpLambda.resize   (NDIM*numPtcls);
  Expm1Lambda.resize (NDIM*numPtcls);
  RArray.resize      (NDIM*numPtcls);
  RpArray.resize     (NDIM*numPtcls);
  TempArray.resize   (NDIM*numPtcls);
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

  DumpRho = in.OpenSection("RhoDump");
  if (DumpRho) {
    int nx, ny, nz;
    assert (in.ReadVar("Nx", nx));
    assert (in.ReadVar("Ny", ny));
    assert (in.ReadVar("Nz", nz));
    Rho.resize(nx, ny, nz);
    in.ReadVar ("DumpBandRho", DumpBandRho, false);
    if (DumpBandRho) 
      if (PathData.Actions.NodalActions(0) != NULL)
	if (PathData.Actions.NodalActions(0)->Type() == GROUND_STATE_FP) {
	  FixedPhaseActionClass &FP = 
	    *((FixedPhaseActionClass*) PathData.Actions.NodalActions(0));
	  int numBands = FP.GetNumBands();
	  BandRho.resize(nx, ny, nz, numBands);
	}
    in.CloseSection (); // "RhoDump"
  }

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

  /// Now initialize Rp
  for (int i=0; i<Rp.size(); i++)
    Rp(i) = R(i) - TimeStep*V(i);
}
