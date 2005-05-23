#include "GroundStateNodalActionClass.h"
#include "../PathDataClass.h"
#include "../Common/MatrixOps/MatrixOps.h"

GroundStateClass::GroundStateClass(PathDataClass &pathData) :
PathData (pathData), Path(pathData.Path)
{
  
}

bool
GroundStateClass::IonsHaveMoved()
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
    cerr << "Ions have moved.\n";
  return changed;
}

double
GroundStateClass::SimpleDistance (int slice, int speciesNum)
{
  double det = GradientDet (slice, speciesNum);
  double gradMag = 0.0;
  for (int i=0; i<Path.Species(speciesNum).NumParticles; i++)
    gradMag += dot (Gradient(i), Gradient(i));

  // Newton Raphson approximation to distance to node 
  double simpleDist = det/sqrt(gradMag);
//   double lineDist = LineSearchDistance(slice, speciesNum);
//   fprintf (stdout, "%1.16e %1.16e\n", simpleDist, lineDist);
  return (simpleDist);
}

double 
GroundStateClass::LineSearchDistance (int slice, int speciesNum)
{
  SpeciesClass &species = Path.Species(speciesNum);
  double epsilon = 1.0e-4 * sqrt (4.0*species.lambda*Path.tau);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;

  //  double maxDist = MaxDist (slice);
  double det0, det;
  int N = last-first+1;
  double retVal;

  // Save the current path
  for (int i=0; i < N; i++)
    Temp(i) = Path(slice,i+first);

  det0 = GradientDet (slice, speciesNum);
  if (det0 < 0.0)
    return -1.0;

  double grad2=0.0;
  for (int i=0; i<N; i++)
    grad2 += dot (Gradient(i), Gradient(i));
  double gradMag = sqrt (grad2);

  for (int i=0; i<N; i++)
    Gradient(i) = (1.0/gradMag)*Gradient(i);

  double dist = det0/gradMag;


  double minFactor, maxFactor, tryFactor, newDet;
  minFactor = 0.0;
  maxFactor = 0.5;

  bool done = false;
  // First, find first sign change
  det = det0;
  double maxDist=sqrt(dot(Path.GetBox(),Path.GetBox()));
  while ((det*det0 > 0.0) && ((maxFactor*dist)<maxDist)) {
    maxFactor *= 2.0;
    for (int i=0; i<N; i++) 
      Path (slice, i+first) = Temp(i) - maxFactor*dist*Gradient(i);
    det = Det(slice, speciesNum);
  }

  if (det*det0 >= 0.0)
    retVal = maxDist;
  else {
    // Now, do a bisection search for the sign change.
    while (((maxFactor-minFactor)*dist > epsilon) 
	   && (minFactor*dist < maxDist)) {
      tryFactor = 0.5*(maxFactor+minFactor);
      for (int i=0; i<N; i++)
	Path (slice, i+first) = Temp(i) - tryFactor*dist*Gradient(i);
      det = Det (slice, speciesNum);
      if (det*det0 > 0.0)
	minFactor = tryFactor;
      else
	maxFactor = tryFactor;
    }
    if (minFactor*dist >= maxDist)
      retVal = maxDist;
    else 
      retVal = dist * tryFactor;
  }

  // Restore original Path position
  for (int i=0; i<N; i++)
    Path (slice, i+first) = Temp(i);
  return retVal; //min (maxDist, retVal);
}




double
GroundStateClass::Action (int slice1, int slice2, 
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

//   for (int i=0; i<activeParticles.size(); i++) 
//     if (Path.ParticleSpeciesNum(activeParticles(i)) == UpSpeciesNum)
//       doUp = true;
//     else if (Path.ParticleSpeciesNum(activeParticles(i)) == DownSpeciesNum)
//       doDown = true;
  doUp   = (speciesNum == UpSpeciesNum);
  doDown = (speciesNum == DownSpeciesNum);
  if (speciesNum == IonSpeciesNum) {
    doUp   = true;
    doDown = true;
  }
  
  if (!(doUp || doDown))
    cerr << "No doing either up or down.  Hmmm...\n";

  double action = 0.0;
  if (doUp) {
    for (int slice=slice1; slice <= slice2; slice++) {
      UpDists(slice) = LineSearchDistance(slice, UpSpeciesNum);
      if (UpDists(slice) <= 0.0) {
	cerr << "negative distance at slice " << slice << endl;
	action += 1.0e100;
      }
    }
    for (int link=slice1; link < slice2; link++) {
      if ((UpDists(link)>0.0) && (UpDists(link+1)>0.0)) {
	double expval = exp(-UpDists(link)*UpDists(link+1)*lambdaTauInv);
	action -= log1p(-expval);
      }
    }
  }
  if (doDown) {
    for (int slice=slice1; slice <= slice2; slice++) {
      DownDists(slice) = LineSearchDistance(slice, DownSpeciesNum);
      if (DownDists(slice) <= 0.0) {
	cerr << "negative distance at slice " << slice << endl;
	action += 1.0e100;
      }
    }
    for (int link=slice1; link < slice2; link++) {
      if ((DownDists(link)>0.0) && (DownDists(link+1)>0.0)) {
	double expval = exp(-DownDists(link)*DownDists(link+1)*lambdaTauInv);
	action -= log1p(-expval);
      }
    }
  }
  return action;
}

double
GroundStateClass::d_dBeta (int slice1, int slice2, int level,
			   int speciesNum)
{
//   bool error = false;
//   for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
//     double upDist = LineSearchDistance (slice, UpSpeciesNum);
//     if (fabs(upDist -UpDists(slice)) > 1.0e-12) {
//       fprintf (stderr, 
// 	       "Cached updist = %1.12e  Actual updist = %1.12e slice=%d\n",
// 	       UpDists(slice), upDist, slice);
//       error = true;
//     }
//     double downDist = LineSearchDistance (slice, DownSpeciesNum);
//     if (fabs(downDist - DownDists(slice))> 1.0e-12){
//       fprintf (stderr, 
// 	       "Cached downdist = %1.12e  Actual downdist = %1.12e slice=%d",
// 	       DownDists(slice), downDist, slice);
//       error = true;
//     }
//   }
//   if (!error)
//     cerr << "No errors detected in GroundStateClass:;d_dBeta.\n";
  double du = 0.0;
  int skip = 1<<level;
  double levelTau =(double)skip * Path.tau;
  double lambdaTauInv = 1.0/(levelTau*Path.Species(UpSpeciesNum).lambda);
  if (speciesNum == UpSpeciesNum)
    for (int slice=slice1; slice<slice2; slice+=skip) {
      double prod = UpDists(slice)*UpDists(slice+skip);
      if (prod <= 0.0)
	cerr << "We have an invalid node-crossing path in "
	     << "GroundStateNodalActionClass::d_dBeta";
      double prod_llt = prod * lambdaTauInv;
      double exp_m1 = expm1 (prod_llt);
      if (isnormal(exp_m1))
	du += prod_llt / (levelTau*exp_m1);
    }
  if (speciesNum == DownSpeciesNum)
    for (int slice=slice1; slice<slice2; slice+=skip) {
      double prod = DownDists(slice)*DownDists(slice+skip);
      if (prod <= 0.0)
	cerr << "We have an invalid node-crossing path in "
	     << "GroundStateNodalActionClass::d_dBeta";
      double prod_llt = prod * lambdaTauInv;
      double exp_m1 = expm1 (prod_llt);
      if (isnormal(exp_m1))
	du += prod_llt / (levelTau*exp_m1);
    }
  return du/(double)Path.TotalNumSlices;
}

void 
GroundStateClass::Read(IOSectionClass &in)
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
  
  Workspace.resize (DetCofactorsWorksize(NumUp));
  Matrix.resize(NumUp, NumUp);
  Cofactors.resize(NumUp, NumUp);
  GradMat.resize(NumUp, NumUp);
  Gradient.resize(NumUp);
  Temp.resize(NumUp);
  UpDists.resize(PathData.NumTimeSlices());
  DownDists.resize(PathData.NumTimeSlices());
  Rions.resize(NumIons);

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
  Array<double,4> initData(nx,ny,nz,NumBands);
  initData = 0.0;
  BandSplines.Init (&xGrid, &yGrid, &zGrid, initData, true);
  //  UpdateBands();
}



double
GroundStateClass::GradientDet(int slice, int speciesNum)
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
  //  cerr << "Analytic gradient = " << Gradient << endl;
//   GradientDetFD(slice, speciesNum);
//   cerr << "FD gradient = " << Gradient << endl;
  return det;
}


double
GroundStateClass::GradientDetFD(int slice, int speciesNum)
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
GroundStateClass::UpdateBands()
{
  // Now, make bands real and put into splines
  Array<double,4> data(xGrid.NumPoints, yGrid.NumPoints, zGrid.NumPoints, NumBands);

  // Only do calcualtion if I'm proc 0
  if (Path.Communicator.MyProc() == 0) {
    cerr << "Updating bands.\n";
    SpeciesClass& ionSpecies = Path.Species(IonSpeciesNum);
    int first = ionSpecies.FirstPtcl;
    for (int i=0; i<NumIons; i++)
      Rions(i) = Path(0,i+first);
    System->SetIons (Rions);
    System->DiagonalizeH();

    for (int band=0; band<NumBands; band++) {
      System->SetRealSpaceBandNum(band);
      complex<double> c0 = 
	System->RealSpaceBand(xGrid.NumPoints/2, 
			      yGrid.NumPoints/2, 
			      zGrid.NumPoints/2);
      double phi = -atan2 (c0.imag(), c0.real());
      complex<double> c(cos(phi), sin(phi));
      for (int ix=0; ix<xGrid.NumPoints-1; ix++)
	for (int iy=0; iy<yGrid.NumPoints-1; iy++)
	  for (int iz=0; iz<zGrid.NumPoints-1; iz++)
	    data(ix,iy,iz,band) = real(c*System->RealSpaceBand(ix,iy,iz));
    }
    MakePeriodic (data);
  }
  Path.Communicator.Broadcast(0, data);
  
  /// DEBUG DEBUG DEBUG DEBUG
  for (int band=0; band<NumBands; band++) {
    char fname[100];
    int proc = 
    snprintf (fname, 100, "band%d-%d.dat", band, Path.Communicator.MyProc());
    FILE *fout = fopen (fname, "w");
    for (int ix=0; ix<xGrid.NumPoints; ix++)
      for (int iy=0; iy<yGrid.NumPoints; iy++)
	for (int iz=0; iz<zGrid.NumPoints; iz++)
	  fprintf (fout, "%1.12e\n", data(ix,iy,iz,band));
    fclose (fout);
  }

  BandSplines.Init (&xGrid, &yGrid, &zGrid, data, true);
}


bool
GroundStateClass::IsPositive (int slice, int speciesNum)
{
  return (Det(slice, speciesNum) > 0.0);  
}

double
GroundStateClass::Det (int slice, int speciesNum)
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
GroundStateClass::GetMatrix (int slice, int speciesNum)
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
GroundStateNodalActionClass::IsPositive (int slice)
{
  return GroundState.IsPositive(slice, SpeciesNum);
}

double
GroundStateNodalActionClass::Det (int slice)
{
  return GroundState.Det(slice, SpeciesNum);
}

Array<double,2>
GroundStateNodalActionClass::GetMatrix (int slice)
{
  return GroundState.GetMatrix(slice, SpeciesNum);
}

double
GroundStateNodalActionClass::Action (int slice1, int slice2, 
				     const Array<int,1> &activeParticles,
				     int level)
{
  double action = GroundState.Action (slice1, slice2, activeParticles, level,
				      SpeciesNum);
//   cerr << "species = " << Path.Species(SpeciesNum).Name 
//        << " NodalAction = " <<  action << endl;
  return action;
}

double
GroundStateNodalActionClass::d_dBeta(int slice1, int slice2, int level)
{
  return GroundState.d_dBeta (slice1, slice2, level, SpeciesNum);
}

void GroundStateClass::ShiftData(int slicesToShift, int speciesNum)
{
  if ((speciesNum == UpSpeciesNum) || (speciesNum==DownSpeciesNum)) {
//     if (speciesNum == UpSpeciesNum) 
//       cerr << "Shifting up species by " << slicesToShift << endl;
//     if (speciesNum == DownSpeciesNum) 
//       cerr << "Shifting down species by " << slicesToShift << endl;
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
GroundStateClass::AcceptCopy (int slice1, int slice2)
{
  UpDists.AcceptCopy (slice1, slice2);
  DownDists.AcceptCopy (slice1, slice2);
}

void 
GroundStateClass::RejectCopy (int slice1, int slice2)
{
  UpDists.RejectCopy (slice1, slice2);
  DownDists.RejectCopy (slice1, slice2);
}

void
GroundStateNodalActionClass::ShiftData (int slices2Shift)
{
  GroundState.ShiftData (slices2Shift, SpeciesNum);
}

void
GroundStateNodalActionClass::AcceptCopy (int slice1, int slice2)
{
  GroundState.AcceptCopy (slice1, slice2);
}

void
GroundStateNodalActionClass::RejectCopy (int slice1, int slice2)
{
  GroundState.RejectCopy (slice1, slice2);
}


void 
GroundStateClass::Init(int speciesNum)
{
  if (speciesNum == UpSpeciesNum) {
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
      UpDists[0](slice) = LineSearchDistance(slice, speciesNum);
      UpDists[1](slice) = UpDists[0](slice);
    }
    cerr << "Initializing up species.\n";
  }

  if (speciesNum == DownSpeciesNum) {
    for (int slice=0; slice<Path.NumTimeSlices(); slice++) {
      DownDists[0](slice) = LineSearchDistance(slice, speciesNum);
      DownDists[1](slice) = DownDists[0](slice);
    }
    cerr << "Initializing down species.\n";
  }
}

void 
GroundStateNodalActionClass::Init()
{
  GroundState.Init (SpeciesNum);
}

bool
GroundStateNodalActionClass::IsGroundState()
{
  return true;
}
