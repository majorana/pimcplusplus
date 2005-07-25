#include "FixedPhaseActionClass.h"
#include "../PathDataClass.h"
#include "../Common/MatrixOps/MatrixOps.h"

FixedPhaseClass::FixedPhaseClass(PathDataClass &pathData) :
PathData (pathData), Path(pathData.Path)
{
  
}

bool
FixedPhaseClass::IonsHaveMoved()
{
  SpeciesClass &species = Path.Species(IonSpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  int ptcl = first;
  bool changed = false;
  for (int ptcl=first; ptcl<=last; ptcl++) {
    dVec diff = Path(0,ptcl) - System->GetIonPos(ptcl-first);
    if (dot (diff,diff) > 1.0e-16)
      changed = true;
  }
  if (changed)
    perr << "Ions have moved.\n";
  return changed;
}




double
FixedPhaseClass::Action (int slice1, int slice2, 
			  const Array<int,1> &activeParticles, 
			  int level, int speciesNum)
{
  double lambdaTauInv = 1.0/(Path.tau*Path.Species(UpSpeciesNum).lambda);

  // The nodal action should only be used at level = 0;
  assert (level == 0);
  bool doUp   = false;
  bool doDown = false;
  if (IonsHaveMoved()) 
    UpdateBands();

  doUp   = (speciesNum == UpSpeciesNum);
  doDown = (speciesNum == DownSpeciesNum);
  if (speciesNum == IonSpeciesNum) {
    doUp   = true;
    doDown = true;
  }
  
  if (!(doUp || doDown))
    perr << "Not doing either up or down.  Hmmm...\n";

  double action = 0.0;
  if (doUp) {
    for (int slice=slice1; slice <= slice2; slice++) 
      CalcVals (slice, UpSpeciesNum);
    for (int link=slice1; link < slice2; link++) 
      action += 0.0;
  }
  if (doDown) {
    for (int slice=slice1; slice <= slice2; slice++) 
      CalcVals (slice, DownSpeciesNum);
    for (int link=slice1; link < slice2; link++) 
      action += 0.0;
  }
  return action;
}

double
FixedPhaseClass::d_dBeta (int slice1, int slice2, int level,
			  int speciesNum)
{

}

void 
FixedPhaseClass::Read(IOSectionClass &in)
{
  string speciesString;
  assert (in.ReadVar ("UpSpecies", speciesString));
  UpSpeciesNum = Path.SpeciesNum (speciesString);
  if (in.ReadVar ("DownSpecies", speciesString))
    DownSpeciesNum = Path.SpeciesNum(speciesString);
  else
    DownSpeciesNum = -1;

  assert (in.ReadVar ("IonSpecies", speciesString));
  IonSpeciesNum = Path.SpeciesNum (speciesString);

  assert (in.ReadVar ("kCut", kCut));
  
  NumIons =  Path.Species(IonSpeciesNum).NumParticles;
  NumUp   =  Path.Species(UpSpeciesNum).NumParticles;
  NumDown = (DownSpeciesNum != -1) ? 
    Path.Species(DownSpeciesNum).NumParticles : -1;

  if (DownSpeciesNum != -1)
    assert (NumUp == NumDown);
  
  Workspace.resize (ComplexDetCofactorsWorksize(NumUp));
  Matrix.resize(NumUp, NumUp);
  Cofactors.resize(NumUp, NumUp);
  GradMat.resize(NumUp, NumUp);
  Gradient.resize(NumUp);
  Temp.resize(NumUp);
  GUp.resize(PathData.NumTimeSlices());
  gUp.resize(PathData.NumTimeSlices());
  vUp.resize(PathData.NumTimeSlices());
  GDown.resize(PathData.NumTimeSlices());
  gDown.resize(PathData.NumTimeSlices());
  vDown.resize(PathData.NumTimeSlices());
  Rions.resize(NumIons);
  Rions = Vec3(0.0, 0.0, 0.0);

  /////////////////////////////////
  // Setup the plane wave system //
  /////////////////////////////////
  NumBands = max(NumUp, NumDown);
  System = new SystemClass (NumBands);
  PH = &PathData.Actions.GetPotential (IonSpeciesNum, UpSpeciesNum);
  Vec3 gamma (0.0, 0.0, 0.0);
  System->Setup (Path.GetBox(), gamma, kCut, *PH);

  /////////////////////////////
  // Setup the ion positions //
  /////////////////////////////
  SpeciesClass& ionSpecies = Path.Species(IonSpeciesNum);
  int first = ionSpecies.FirstPtcl;
  for (int i=0; i<NumIons; i++)
    Rions(i) = Path(0,i+first);
  Rions = Vec3(0.0, 0.0, 0.0);
  System->SetIons (Rions);

  ////////////////////////////////////////
  // Setup real space grids and splines //
  ////////////////////////////////////////
  int nx, ny, nz;
  System->GetBoxDims (nx, ny, nz);
  // Increment sizes to compensate for PBC
  nx++; ny++; nz++;
  Vec3 box;
  box = Path.GetBox();
  xGrid.Init (-0.5*box[0], 0.5*box[0], nx);
  yGrid.Init (-0.5*box[1], 0.5*box[1], ny);
  zGrid.Init (-0.5*box[2], 0.5*box[2], nz);
  Array<complex<double>,4> initData(nx,ny,nz,NumBands);
  initData = 0.0;
  BandSplines.Init (&xGrid, &yGrid, &zGrid, initData, true);
  //  UpdateBands();
}



double
FixedPhaseClass::GradientDet(int slice, int speciesNum)
{
  SpeciesClass &species = Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  // First, fill up determinant matrix
  Array<double,1> vals;
  Array<Vec3,1> grads;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    vals.reference(Matrix(j,Range::all()));
    grads.reference(GradMat(j,Range::all()));
    BandSplines.ValGrad(r_j[0], r_j[1], r_j[2], vals, grads);
  }

  // Compute determinant and cofactors
  Cofactors = Matrix;
  double det = DetCofactors (Cofactors, Workspace);

  // Now, compute gradient
  Gradient = Vec3(0.0, 0.0, 0.0);
  for (int i=0; i<N; i++) {
//     Vec3 r_i = Path(slice, first+i);
//     BandSplines.Grad(r_i[0], r_i[1], r_i[2], Temp);
    for (int j=0; j<N; j++)
      Gradient(i) += Cofactors(i,j)*GradMat(i,j);
  }
  //  perr << "Analytic gradient = " << Gradient << endl;
//   GradientDetFD(slice, speciesNum);
//   perr << "FD gradient = " << Gradient << endl;
  return det;
}


double
FixedPhaseClass::GradientDetFD(int slice, int speciesNum)
{
  SpeciesClass &species = Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  const double eps = 1.0e-6;

  // First, fill up determinant matrix
  Array<double,1> row;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    row.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], row);
  }

  double det = Determinant (Matrix);

  for (int i=0; i<N; i++) {
    Vec3 plus, minus;
    plus = 0.0; minus = 0.0;
    Vec3 r_i = Path(slice, first+i);
    Vec3 r;
    row.reference(Matrix(i,Range::all()));
    for (int dim=0; dim<NDIM; dim++) {
      r = r_i;
      r[dim] += eps;
      Path.PutInBox(r);
      BandSplines(r[0], r[1], r[2], row);
      plus[dim] = Determinant (Matrix);
      
      r = r_i;
      r[dim] -= eps;
      Path.PutInBox(r);
      BandSplines(r[0], r[1], r[2], row);
      minus[dim] = Determinant (Matrix);
    }
    // Reset row to orignal value
    r = r_i;
    Path.PutInBox(r);
    BandSplines(r[0], r[1], r[2], row);
    // And compute gradient
    Gradient(i) = (plus-minus)/(2.0*eps);
  }
  return det;
}


void
FixedPhaseClass::UpdateBands()
{
  // Now, make bands real and put into splines
  Array<double,4> data(xGrid.NumPoints, yGrid.NumPoints, zGrid.NumPoints, NumBands);
  SpeciesClass& ionSpecies = Path.Species(IonSpeciesNum);
  int first = ionSpecies.FirstPtcl;
  for (int i=0; i<NumIons; i++)
    Rions(i) = Path(0,i+first);
  System->SetIons (Rions);

  // Only do calcualtion if I'm proc 0
  if (Path.Communicator.MyProc() == 0) {
    perr << "Updating bands.\n";
    System->DiagonalizeH();

    for (int band=0; band<NumBands; band++) {
      System->SetRealSpaceBandNum(band);
      for (int ix=0; ix<xGrid.NumPoints-1; ix++)
	for (int iy=0; iy<yGrid.NumPoints-1; iy++)
	  for (int iz=0; iz<zGrid.NumPoints-1; iz++)
	    data(ix,iy,iz,band) = System->RealSpaceBand(ix,iy,iz);
    }
    MakePeriodic (data);
  }
  Path.Communicator.Broadcast(0, data);
  
  /// DEBUG DEBUG DEBUG DEBUG
//   for (int band=0; band<NumBands; band++) {
//     char fname[100];
//     int proc = 
//     snprintf (fname, 100, "band%d-%d.dat", band, Path.Communicator.MyProc());
//     FILE *fout = fopen (fname, "w");
//     for (int ix=0; ix<xGrid.NumPoints; ix++)
//       for (int iy=0; iy<yGrid.NumPoints; iy++)
// 	for (int iz=0; iz<zGrid.NumPoints; iz++)
// 	  fprintf (fout, "%1.12e\n", data(ix,iy,iz,band));
//     fclose (fout);
//   }

  BandSplines.Init (&xGrid, &yGrid, &zGrid, data, true);
}


bool
FixedPhaseClass::IsPositive (int slice, int speciesNum)
{
  return (Det(slice, speciesNum) > 0.0);  
}

double
FixedPhaseClass::Det (int slice, int speciesNum)
{
  if (IonsHaveMoved()) 
    UpdateBands();

  SpeciesClass& species = Path.Species(speciesNum);
  int first = species.FirstPtcl;
  int N = species.NumParticles;
  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }
  return Determinant (Matrix);
}


Array<double,2>
FixedPhaseClass::GetMatrix (int slice, int speciesNum)
{
  if (IonsHaveMoved()) 
    UpdateBands();

  SpeciesClass& species = Path.Species(speciesNum);
  int first = species.FirstPtcl;
  int N = species.NumParticles;
  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }
  return Matrix;
}


bool
FixedPhaseActionClass::IsPositive (int slice)
{
  return FixedPhase.IsPositive(slice, SpeciesNum);
}

double
FixedPhaseActionClass::Det (int slice)
{
  return FixedPhase.Det(slice, SpeciesNum);
}

Array<double,2>
FixedPhaseActionClass::GetMatrix (int slice)
{
  return FixedPhase.GetMatrix(slice, SpeciesNum);
}

double
FixedPhaseActionClass::Action (int slice1, int slice2, 
				     const Array<int,1> &activeParticles,
				     int level)
{
  double action = FixedPhase.Action (slice1, slice2, activeParticles, level,
				      SpeciesNum);
  return action;
}

double
FixedPhaseActionClass::d_dBeta(int slice1, int slice2, int level)
{
  return FixedPhase.d_dBeta (slice1, slice2, level, SpeciesNum);
}

void FixedPhaseClass::ShiftData(int slicesToShift, int speciesNum)
{
  if ((speciesNum == UpSpeciesNum) || (speciesNum==DownSpeciesNum)) {
//     if (speciesNum == UpSpeciesNum) 
//       perr << "Shifting up species by " << slicesToShift << endl;
//     if (speciesNum == DownSpeciesNum) 
//       perr << "Shifting down species by " << slicesToShift << endl;
    Mirrored1DClass<double>& dists=
      ((speciesNum==UpSpeciesNum) ? UpDists : DownDists);
    CommunicatorClass &comm = Path.Communicator;

    int numProcs=comm.NumProcs();
    int myProc=comm.MyProc();
    int recvProc, sendProc;
    int numSlices  = Path.NumTimeSlices();
    assert(abs(slicesToShift)<numSlices);
    sendProc=(myProc+1) % numProcs;
    recvProc=((myProc-1) + numProcs) % numProcs;
    if (slicesToShift<0)
      swap (sendProc, recvProc);
    
    /// First shifts the data in the A copy left or right by the
    /// appropriate amount 
    if (slicesToShift>0)
      for (int slice=numSlices-1; slice>=slicesToShift;slice--)
	dists[0](slice)=dists[0](slice-slicesToShift);
    else 
      for (int slice=0; slice<numSlices+slicesToShift;slice++)
	dists[0](slice)=dists[0](slice-slicesToShift);
    
    
    /// Now bundle up the data to send to adjacent processor
    int bufferSize=abs(slicesToShift);
    Array<double,1> sendBuffer(bufferSize), receiveBuffer(bufferSize);
    int startSlice;
    int buffIndex=0;
    if (slicesToShift>0) {
      startSlice=numSlices-slicesToShift;
      for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++) {
	/// If shifting forward, don't send the last time slice (so always)
	/// send slice-1
	sendBuffer(buffIndex)=dists[1](slice-1);
	buffIndex++;
      }
    }
    else {
      startSlice=0;
      for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
	/// If shifting backward, don't send the first time slice (so always)
	/// send slice+1
	sendBuffer(buffIndex)=dists[1](slice+1);
	buffIndex++;
      }
    }
    
    /// Send and receive data to/from neighbors.
    comm.SendReceive(sendProc, sendBuffer,recvProc, receiveBuffer);
    
    if (slicesToShift>0)
      startSlice=0;
    else 
      startSlice=numSlices+slicesToShift;
    
    /// Copy the data into the A copy
    buffIndex=0;
    for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
      dists[0](slice)=receiveBuffer(buffIndex);
      buffIndex++;
    }
    
    // Now copy A into B, since A has all the good, shifted data now.
    for (int slice=0; slice<numSlices; slice++)
      dists[1](slice) = dists[0](slice);
    // And we're done! 
  } 
}

void 
FixedPhaseClass::AcceptCopy (int slice1, int slice2)
{
  UpDists.AcceptCopy (slice1, slice2);
  DownDists.AcceptCopy (slice1, slice2);
}

void 
FixedPhaseClass::RejectCopy (int slice1, int slice2)
{
  UpDists.RejectCopy (slice1, slice2);
  DownDists.RejectCopy (slice1, slice2);
}

void
FixedPhaseActionClass::ShiftData (int slices2Shift)
{
  FixedPhase.ShiftData (slices2Shift, SpeciesNum);
}

void
FixedPhaseActionClass::AcceptCopy (int slice1, int slice2)
{
  FixedPhase.AcceptCopy (slice1, slice2);
}

void
FixedPhaseActionClass::RejectCopy (int slice1, int slice2)
{
  FixedPhase.RejectCopy (slice1, slice2);
}


void 
FixedPhaseClass::Init(int speciesNum)
{
  if (speciesNum == UpSpeciesNum) {
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
      UpDists[0](slice) = LineSearchDistance(slice, speciesNum);
      UpDists[1](slice) = UpDists[0](slice);
    }
    perr << "Initializing up species.\n";
  }

  if (speciesNum == DownSpeciesNum) {
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
      DownDists[0](slice) = LineSearchDistance(slice, speciesNum);
      DownDists[1](slice) = DownDists[0](slice);
    }
    perr << "Initializing down species.\n";
  }
}

void 
FixedPhaseActionClass::Init()
{
  FixedPhase.Init (SpeciesNum);
}

bool
FixedPhaseActionClass::IsGroundState()
{
  return true;
}
