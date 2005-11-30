#include "FixedPhaseActionClass.h"
#include "../PathDataClass.h"
#include <Common/MatrixOps/MatrixOps.h>

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

inline Vec3 real(cVec3 v)
{
  return Vec3(v[0].real(), v[1].real(), v[2].real());
}

inline Vec3 imag(cVec3 v)
{
  return Vec3(v[0].imag(), v[1].imag(), v[2].imag());
}

inline double mag2 (complex<double> z)
{
  return (z.real()*z.real()+z.imag()*z.imag());
}

double
FixedPhaseClass::CalcGrad2 (int slice, int species)
{
  int N = Path.Species(species).NumParticles;
  /// calculate \f$\det|u|\f$ and \f$ \nabla det|u| \f$.
  complex<double> detu = GradientDet(slice, species);
  double detu2 = mag2(detu);
  double detu2Inv = 1.0/detu2;
  
  double grad2 = 0.0;
  for (int i=0; i<N; i++) {
    Vec3 grad = 
      (detu.real()*imag(Gradient(i)) - detu.imag()*real(Gradient(i)))*detu2Inv - kVec;
    grad2 += dot(grad,grad);
  }
  return grad2;
}

/// This function returns the calculates the gradient squared for each
/// particle in the species.
void
FixedPhaseClass::CalcGrad2 (int slice, int species, 
			    Array<double,1> &grad2)
{
  int N = Path.Species(species).NumParticles;
  if (grad2.size() != N)
    grad2.resize(N);
  /// calculate \f$\det|u|\f$ and \f$ \nabla det|u| \f$.
  complex<double> detu = GradientDet(slice, species);
  double detu2 = mag2(detu);
  double detu2Inv = 1.0/detu2;
  
  for (int i=0; i<N; i++) {
    Vec3 grad = (detu.real()*imag(Gradient(i)) - 
		 detu.imag()*real(Gradient(i)))*detu2Inv - kVec;
    grad2(i) = dot(grad,grad);
  }
}


double
FixedPhaseClass::Action (int slice1, int slice2, 
			 const Array<int,1> &activeParticles, 
			 int level, int speciesNum)
{
  int skip = 1<<level;
  double levelTau = ldexp(Path.tau, level);
  const double lambda = Path.Species(UpSpeciesNum).lambda;

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
    perr << "Not doing either up or down.  Hmmm... " 
	 << " speciesNum = " << speciesNum << endl;

  double action = 0.0;
  if (doUp) {
    for (int slice=slice1; slice <= slice2; slice+=skip) 
      UpGrad2(slice) = CalcGrad2 (slice, UpSpeciesNum);
    for (int link=slice1; link < slice2; link+=skip) 
      action += 0.5*lambda*levelTau*(UpGrad2(link)+UpGrad2(link+skip));
  }
  if (doDown) {
    for (int slice=slice1; slice <= slice2; slice+=skip) 
      DownGrad2(slice) = CalcGrad2 (slice, DownSpeciesNum);
    for (int link=slice1; link < slice2; link+=skip) 
      action += 0.5*lambda*levelTau*(DownGrad2(link)+DownGrad2(link+skip));
  }
  return action;
}

double
FixedPhaseClass::d_dBeta (int slice1, int slice2, int level,
			  int speciesNum)
{
  int skip = 1<<level;
  const double lambda = Path.Species(UpSpeciesNum).lambda;

  // The nodal action should only be used at level = 0;
  assert (level == 0);
  bool doUp   = false;
  bool doDown = false;
  if (IonsHaveMoved()) 
    UpdateBands();

  doUp   = (speciesNum == UpSpeciesNum);
  doDown = (speciesNum == DownSpeciesNum);
  // I think the follow was causing the nodal energy to be
  // double-counted:
  //   if (speciesNum == IonSpeciesNum) {
  //     doUp   = true;
  //     doDown = true;
  //   }

  double dU = 0.0;
  if (doUp) 
    for (int link=slice1; link < slice2; link+=skip) 
      dU += 0.5*lambda*(UpGrad2(link)+UpGrad2(link+skip));
  if (doDown) 
    for (int link=slice1; link < slice2; link+=skip)
      dU += 0.5*lambda*(DownGrad2(link)+DownGrad2(link+skip));
  return dU;
}

void 
FixedPhaseClass::Read(IOSectionClass &in)
{
#if NDIM==3
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
  /// Note that the twist angles are defined between 0 and 1.  Hence,
  /// k_x = pi*twist_angle_x/L_x
  Array<double,2> twistAngles;
  if (in.ReadVar ("TwistAngles", twistAngles)) {
    if (PathData.IntraComm.MyProc()==0)
      if (twistAngles.extent(0) != PathData.InterComm.NumProcs()) {
	cerr << "Error:  Number of twist angles must match number of clones.\n" 
	     << "        Number of twist angles:  "  << twistAngles.size() << endl 
	     << "        Number of clones:        "  << PathData.InterComm.NumProcs() << endl;
	abort();
      }
    assert (twistAngles.extent(1) == 3);
    kVec[0] = twistAngles(PathData.GetCloneNum(),0) * M_PI/PathData.Path.GetBox()[0];
    kVec[1] = twistAngles(PathData.GetCloneNum(),1) * M_PI/PathData.Path.GetBox()[1];
    kVec[2] = twistAngles(PathData.GetCloneNum(),2) * M_PI/PathData.Path.GetBox()[2];
  }
  else {
    Array<double,1> tmpkVec;
    assert(in.ReadVar ("kVec", tmpkVec));
    kVec = Vec3 (tmpkVec(0), tmpkVec(1), tmpkVec(2));
  }
  
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
  UpGrad2.resize(PathData.NumTimeSlices());
  DownGrad2.resize(PathData.NumTimeSlices());
  Rions.resize(NumIons);
  Rions = Vec3(0.0, 0.0, 0.0);

  /////////////////////////////////
  // Setup the plane wave system //
  /////////////////////////////////
  NumBands = max(NumUp, NumDown);
  System = new MPISystemClass (NumBands, Path.Communicator);
  PH = &PathData.Actions.GetPotential (IonSpeciesNum, UpSpeciesNum);
  //  Vec3 gamma (0.0, 0.0, 0.0);
  System->Setup (Path.GetBox(), kVec, kCut, *PH);

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
#endif
}



complex<double>
FixedPhaseClass::GradientDet(int slice, int speciesNum)
{
#if NDIM==3
  SpeciesClass &species = Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  if (N != Matrix.rows()) {
    cerr << "N = " << N << endl;
    cerr << "Matrix.rows() = " << Matrix.rows() << endl;
  }
  assert (N == Matrix.rows());

  // First, fill up determinant matrix
  Array<complex<double>,1> vals;
  Array<cVec3,1> grads;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    vals.reference(Matrix(j,Range::all()));
    grads.reference(GradMat(j,Range::all()));
    BandSplines.ValGrad(r_j[0], r_j[1], r_j[2], vals, grads);
//     // New way of adding twist term
//     double phi = -dot (kVec, r_j);
//     double sinphi, cosphi;
//     sincos (phi, &sinphi, &cosphi);
//     complex<double> e2iphi = complex<double>(cosphi, sinphi);
//     cVec3 ikVec;
//     ikVec[0] = complex<double>(0.0, kVec[0]);
//     ikVec[1] = complex<double>(0.0, kVec[1]);
//     ikVec[2] = complex<double>(0.0, kVec[2]);
//     for (int k=0; k<N; k++) {
//       grads(k) = e2iphi*(grads(k) - vals(k)* ikVec);
//       vals(k) *= e2iphi;
//    }
  }

  // Compute determinant and cofactors
  Cofactors = Matrix;
  complex<double> det = ComplexDetCofactors (Cofactors, Workspace);

  // Now, compute gradient
  Gradient = cVec3(0.0, 0.0, 0.0);
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++)
      Gradient(i) += Cofactors(i,j)*GradMat(i,j);
  }
  //  perr << "Analytic gradient = " << Gradient << endl;
//   GradientDetFD(slice, speciesNum);
//   perr << "FD gradient = " << Gradient << endl;
  return det;
#endif
}


complex<double>
FixedPhaseClass::GradientDetFD(int slice, int speciesNum)
{
#if NDIM==3
  SpeciesClass &species = Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  const double eps = 1.0e-6;

  // First, fill up determinant matrix
  Array<complex<double>,1> row;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    row.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], row);
  }

  complex<double> det = Determinant (Matrix);

  for (int i=0; i<N; i++) {
    cVec3 plus, minus;
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
#endif
}


void
FixedPhaseClass::UpdateBands()
{
#if NDIM==3
  // Now, make bands real and put into splines
  Array<complex<double>,4> data(xGrid.NumPoints, 
				yGrid.NumPoints, 
				zGrid.NumPoints, NumBands);
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
#endif
}


bool
FixedPhaseClass::IsPositive (int slice, int speciesNum)
{
  if (IonsHaveMoved()) 
    UpdateBands();


  // HACK HACK HACK
//   const int N=200;
//   Array<double,3> grad2Array(N, N, N);
//   Vec3 pos;
//   for (int ix=0; ix<N; ix++) {
//     pos[0] = ((double)ix/(N-1) -0.5)*Path.GetBox()[0];
//     for (int iy=0; iy<N; iy++){
//       pos[1] = ((double)iy/(N-1) -0.5)*Path.GetBox()[1];
//       for (int iz=0; iz<N; iz++) {
// 	pos[2] = ((double)iz/(N-1) -0.5)*Path.GetBox()[2];
// 	Path(slice,Path.Species(speciesNum).FirstPtcl) = pos;
// 	grad2Array(ix,iy,iz) = CalcGrad2 (slice, speciesNum);
//       }
//     }
//   }
//   IOSectionClass out;
//   out.NewFile ("PhaseDump.h5");
//   out.WriteVar("Grad2", grad2Array);
//   out.CloseFile();
   
  double grad2 = CalcGrad2 (slice, speciesNum);
  //  cerr << "grad2 = " << grad2 << endl;
  double lambda = Path.Species(speciesNum).lambda;

  /// HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK HACK
  //return ((lambda*grad2) < 100.0);
  return true;
}

complex<double>
FixedPhaseClass::Det (int slice, int speciesNum)
{
#if NDIM==3
  if (IonsHaveMoved()) 
    UpdateBands();

  SpeciesClass& species = Path.Species(speciesNum);
  int first = species.FirstPtcl;
  int N = species.NumParticles;
  // First, fill up determinant matrix
  Array<complex<double>,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    Path.PutInBox(r_j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }
  return Determinant (Matrix);
#endif
}


// Array<double,2>
// FixedPhaseClass::GetMatrix (int slice, int speciesNum)
// {
//   if (IonsHaveMoved()) 
//     UpdateBands();

//   SpeciesClass& species = Path.Species(speciesNum);
//   int first = species.FirstPtcl;
//   int N = species.NumParticles;
//   // First, fill up determinant matrix
//   Array<double,1> vals;
//   for (int j=0; j<N; j++) {
//     Vec3 r_j = Path(slice, first+j);
//     Path.PutInBox(r_j);
//     vals.reference(Matrix(j,Range::all()));
//     BandSplines(r_j[0], r_j[1], r_j[2], vals);
//   }
//   return Matrix;
// }


bool
FixedPhaseActionClass::IsPositive (int slice)
{
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.IsPositive(slice, SpeciesNum);
  else
    return FixedPhaseB.IsPositive(slice, SpeciesNum);
}

complex<double>
FixedPhaseActionClass::Det (int slice)
{
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.Det(slice, SpeciesNum);
  else
    return FixedPhaseB.Det(slice, SpeciesNum);
}

double
FixedPhaseActionClass::SingleAction (int slice1, int slice2, 
				     const Array<int,1> &activeParticles,
				     int level)
{
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.Action (slice1, slice2, activeParticles, level, SpeciesNum);
  else
    return FixedPhaseB.Action (slice1, slice2, activeParticles, level, SpeciesNum);
}

double
FixedPhaseActionClass::d_dBeta(int slice1, int slice2, int level)
{
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.d_dBeta (slice1, slice2, level, SpeciesNum);
  else
    return FixedPhaseB.d_dBeta (slice1, slice2, level, SpeciesNum);
}

void FixedPhaseClass::ShiftData(int slicesToShift, int speciesNum)
{
  if ((speciesNum == UpSpeciesNum) || (speciesNum==DownSpeciesNum)) {
    Mirrored1DClass<double>& grad2 =
      ((speciesNum==UpSpeciesNum) ? UpGrad2 : DownGrad2);
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
	grad2[0](slice)=grad2[0](slice-slicesToShift);
    else 
      for (int slice=0; slice<numSlices+slicesToShift;slice++)
	grad2[0](slice)=grad2[0](slice-slicesToShift);
    
    
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
	sendBuffer(buffIndex)=grad2[1](slice-1);
	buffIndex++;
      }
    }
    else {
      startSlice=0;
      for (int slice=startSlice; slice<startSlice+abs(slicesToShift);slice++){
	/// If shifting backward, don't send the first time slice (so always)
	/// send slice+1
	sendBuffer(buffIndex)=grad2[1](slice+1);
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
      grad2[0](slice)=receiveBuffer(buffIndex);
      buffIndex++;
    }
    
    // Now copy A into B, since A has all the good, shifted data now.
    for (int slice=0; slice<numSlices; slice++)
      grad2[1](slice) = grad2[0](slice);
    // And we're done! 
  } 
}

void 
FixedPhaseClass::AcceptCopy (int slice1, int slice2)
{
  UpGrad2.AcceptCopy(slice1, slice2);
  DownGrad2.AcceptCopy(slice1, slice2);
}

void 
FixedPhaseClass::RejectCopy (int slice1, int slice2)
{
  UpGrad2.RejectCopy(slice1, slice2); 
  DownGrad2.RejectCopy(slice1, slice2); 
}

void
FixedPhaseActionClass::ShiftData (int slices2Shift)
{
  FixedPhaseA.ShiftData (slices2Shift, SpeciesNum);
  FixedPhaseB.ShiftData (slices2Shift, SpeciesNum);
}

void
FixedPhaseActionClass::AcceptCopy (int slice1, int slice2)
{
  FixedPhaseA.AcceptCopy (slice1, slice2);
  FixedPhaseB.AcceptCopy (slice1, slice2);
}

void
FixedPhaseActionClass::RejectCopy (int slice1, int slice2)
{
  FixedPhaseA.RejectCopy (slice1, slice2);
  FixedPhaseB.RejectCopy (slice1, slice2);
}




void 
FixedPhaseClass::Init(int speciesNum)
{
  SetMode(NEWMODE);

  if (speciesNum == UpSpeciesNum) {
    perr << "Initializing up species.\n";  
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) 
      CalcGrad2(slice, UpSpeciesNum);
    AcceptCopy (0, Path.NumTimeSlices()-1);
  }
  else if (speciesNum == DownSpeciesNum) {
    perr << "Initializing down species.\n";
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) 
      CalcGrad2(slice, DownSpeciesNum);
    AcceptCopy (0, Path.NumTimeSlices()-1);
  }
}

void 
FixedPhaseActionClass::Init()
{
  FixedPhaseA.Init (SpeciesNum);
  FixedPhaseB.Init (SpeciesNum);
}

bool
FixedPhaseActionClass::IsGroundState()
{
  return true;
}


NodeType
FixedPhaseActionClass::Type()
{
  return GROUND_STATE_FP;
}

void
FixedPhaseActionClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "GROUND_STATE_FP");
  Array<double,1> outkVec(3);
  outkVec(0) = Getk()[0];  outkVec(1) = Getk()[1];  outkVec(2) = Getk()[2];
  out.WriteVar ("kVec", outkVec);
  double energy = 0.0;
  if ((SpeciesNum==FixedPhaseA.UpSpeciesNum) || (SpeciesNum == FixedPhaseA.DownSpeciesNum))
    for (int i=0; i<FixedPhaseA.System->GetNumBands(); i++) {
      cerr << "Band Energy = " << FixedPhaseA.System->GetEnergy(i) << endl;
      energy += FixedPhaseA.System->GetEnergy(i);
    }
  out.WriteVar ("Energy", energy);
//   int nx = FixedPhase.BandSplines.Nx;
//   int ny = FixedPhase.BandSplines.Ny;
//   int nz = FixedPhase.BandSplines.Nz;
//   int n  = FixedPhase.BandSplines.N;
//   Array<double,4> Bands(nx, ny, nz, 2*n);
//   for (int ix=0; ix<nx; ix++)
//     for (int iy=0; iy<ny; iy++)
//       for (int iz=0; iz<nz; iz++)
// 	for (int i=0; i<n; i++) {
// 	  Bands(ix,iy,iz,2*i)   = FixedPhase.BandSplines(ix,iy,iz,i).real();
// 	  Bands(ix,iy,iz,2*i+1) = FixedPhase.BandSplines(ix,iy,iz,i).imag();
// 	}
//   out.WriteVar("Bands", Bands);
//   Vec3 r;
//   for (int ix=0; ix<nx; ix++) {
//     r[0] = (*FixedPhase.BandSplines.Xgrid)(ix);
//     for (int iy=0; iy<ny; iy++) {
//       r[1] = (*FixedPhase.BandSplines.Ygrid)(iy);
//       for (int iz=0; iz<nz; iz++) {
// 	r[2] = (*FixedPhase.BandSplines.Zgrid)(iz);
// 	double phi = dot(Getk(), r);
// 	for (int i=0; i<n; i++) {
// 	  Bands(ix,iy,iz,2*i)   = cos(phi);
// 	  Bands(ix,iy,iz,2*i+1) = sin(phi);
// 	}
//       }
//     }
//   }
//   out.WriteVar("e2ikr", Bands);
  out.FlushFile();
}

void
FixedPhaseActionClass::Setk(Vec3 k)
{
  FixedPhaseA.Setk(k);
  FixedPhaseB.Setk(k);
}


void
FixedPhaseActionClass::CalcGrad2 (int slice, Array<double,1> &grad2)
{
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.CalcGrad2(slice, SpeciesNum, grad2); 
  else
    return FixedPhaseB.CalcGrad2(slice, SpeciesNum, grad2); 
}

double 
FixedPhaseActionClass::CalcGrad2 (int slice) 
{ 
  if (PathData.Path.GetConfig() == 0)
    return FixedPhaseA.CalcGrad2(slice, SpeciesNum); 
  else
    return FixedPhaseB.CalcGrad2(slice, SpeciesNum); 
}


void
FixedPhaseClass::Setk(Vec3 k)
{
  kVec = k;
  System->Setk(k);
  UpdateBands();
}


