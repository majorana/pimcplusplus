#include "GroundStateNodalActionClass.h"
#include "../PathDataClass.h"
#include "../Common/MatrixOps/MatrixOps.h"

bool
GroundStateNodalActionClass::IonsHaveMoved()
{
  SpeciesClass &species = PathData.Path.Species(IonSpeciesNum);
  int first = species.FirstPtcl;
  int last = species.LastPtcl;
  int ptcl = first;
  bool changed = false;
  for (int ptcl=first; ptcl<=last; ptcl++) {
    if (!(PathData.Path(0,ptcl) == System->GetIonPos(ptcl-first)))
      changed = true;
  }
  return changed;
}

double
GroundStateNodalActionClass::SimpleDistance (int slice, int speciesNum)
{
  double det = GradientDet (slice, speciesNum);
  double gradMag = 0.0;
  for (int i=0; i<Path.Species(speciesNum).NumParticles; i++)
    gradMag += dot (Gradient(i), Gradient(i));

  // Newton Raphson approximation to distance to node 
  return (det/sqrt(gradMag));
}

double 
GroundStateNodalActionClass::LineSearchDistance (int slice, int speciesNum)
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

  GradientDet (slice, speciesNum);
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
  double maxDist=sqrt(dot(PathData.Path.GetBox(),PathData.Path.GetBox()));
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
GroundStateNodalActionClass::Action (int slice1, int slice2, 
				     const Array<int,1> &activeParticles, 
				     int level)
{
  double lambdaTauInv = 1.0/(Path.tau*Path.Species(UpSpeciesNum).lambda);

  // The nodal action should only be used at level = 0;
  assert (level == 0);
  bool doUp   = false;
  bool doDown = false;
  if (IonsHaveMoved()) 
    UpdateBands();

  for (int i=0; i<activeParticles.size(); i++) 
    if (Path.ParticleSpeciesNum(activeParticles(i)) == UpSpeciesNum)
      doUp = true;
    else if (Path.ParticleSpeciesNum(activeParticles(i)) == DownSpeciesNum)
      doDown = true;

  double action = 0.0;
  if (doUp) {
    for (int slice=slice1; slice <= slice2; slice++) {
      UpDists(slice) = SimpleDistance(slice, UpSpeciesNum);
      if (UpDists(slice) < 0.0)
	return 1.0e100;
    }
    for (int link=slice1; link < slice2; link++) 
      action -= log1p(-exp(-UpDists(link)*UpDists(link+1)*lambdaTauInv));
  }
  if (doDown) {
    for (int slice=slice1; slice <= slice2; slice++) {
      double det = GradientDet (slice, UpSpeciesNum);
      double gradMag = 0.0;
      for (int i=0; i<Path.Species(UpSpeciesNum).NumParticles; i++)
	gradMag += dot (Gradient(i), Gradient(i));
      // Newton Raphson approximation to distance to node
      DownDists(slice) = det/sqrt(gradMag);
      if (DownDists(slice) < 0.0)
	return 1.0e100;
    }
    for (int link=slice1; link < slice2; link++) 
      action -= log1p(-exp(-DownDists(link)*DownDists(link+1)*lambdaTauInv));
  }
  return action;
}

double
GroundStateNodalActionClass::d_dBeta (int slice1, int slice2, int level)
{
  return 0.0;
}

void 
GroundStateNodalActionClass::Read(IOSectionClass &in)
{
  string speciesString;
  assert (in.ReadVar ("UpSpecies", speciesString));
  UpSpeciesNum = PathData.Path.SpeciesNum (speciesString);
  if (in.ReadVar ("DownSpecies", speciesString))
    DownSpeciesNum = PathData.Path.SpeciesNum(speciesString);
  else
    DownSpeciesNum = -1;

  assert (in.ReadVar ("IonSpecies", speciesString));
  IonSpeciesNum = PathData.Path.SpeciesNum (speciesString);

  assert (in.ReadVar ("kCut", kCut));
  
  NumIons =  PathData.Path.Species(IonSpeciesNum).NumParticles;
  NumUp   =  PathData.Path.Species(UpSpeciesNum).NumParticles;
  NumDown = (DownSpeciesNum != -1) ? 
    PathData.Path.Species(DownSpeciesNum).NumParticles : -1;

  if (DownSpeciesNum != -1)
    assert (NumUp == NumDown);
  
  Workspace.resize (DetCofactorsWorksize(NumUp));
  Matrix.resize(NumUp, NumUp);
  Cofactors.resize(NumUp, NumUp);
  Gradient.resize(NumUp);
  Temp.resize(NumUp);
  UpDists.resize(PathData.NumTimeSlices());
  DownDists.resize(PathData.NumTimeSlices());
  
  /////////////////////////////////
  // Setup the plane wave system //
  /////////////////////////////////
  NumBands = max(NumUp, NumDown);
  System = new SystemClass (NumBands);
  PH = &PathData.Actions.GetPotential (IonSpeciesNum, UpSpeciesNum);
  Vec3 gamma (0.0, 0.0, 0.0);
  System->Setup (PathData.Path.GetBox(), gamma, kCut, *PH);

  ////////////////////////////////////////
  // Setup real space grids and splines //
  ////////////////////////////////////////
  int nx, ny, nz;
  System->GetBoxDims (nx, ny, nz);
  // Increment sizes to compensate for PBC
  nx++; ny++; nz++;
  Vec3 box;
  box = PathData.Path.GetBox();
  xGrid.Init (-0.5*box[0], 0.5*box[0], nx);
  yGrid.Init (-0.5*box[1], 0.5*box[1], nx);
  zGrid.Init (-0.5*box[2], 0.5*box[2], nx);
  Array<double,4> initData(nx,ny,nz,NumBands);
  initData = 0.0;
  BandSplines.Init (&xGrid, &yGrid, &zGrid, initData, true);
}

double
GroundStateNodalActionClass::Det(int slice, int speciesNum)
{
  SpeciesClass &species = PathData.Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = PathData.Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }

  // Compute determinant and cofactors
  Cofactors = Matrix;
  double det = Determinant (Cofactors);

  return det;
}


double
GroundStateNodalActionClass::GradientDet(int slice, int speciesNum)
{
  SpeciesClass &species = PathData.Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = PathData.Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }

  // Compute determinant and cofactors
  Cofactors = Matrix;
  double det = DetCofactors (Cofactors, Workspace);

  // Now, compute gradient
  Gradient = Vec3(0.0, 0.0, 0.0);
  for (int i=0; i<N; i++) {
    Vec3 r_i = PathData.Path(slice, first+i);
    BandSplines.Grad(r_i[0], r_i[1], r_i[2], Temp);
    for (int j=0; j<N; j++)
      Gradient(j) += Cofactors(i,j)*Temp(j);
  }
  cerr << "Analytic gradient = " << Gradient << endl;
  GradientDetFD(slice, speciesNum);
  cerr << "FD gradient = " << Gradient << endl;
  return det;
}


double
GroundStateNodalActionClass::GradientDetFD(int slice, int speciesNum)
{
  SpeciesClass &species = PathData.Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  const double eps = 1.0e-6;

  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = PathData.Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }

  double det = Determinant (Matrix);

  for (int i=0; i<N; i++) {
    Vec3 plus, minus;
    plus = 0.0; minus = 0.0;
    for (int dim=0; dim<NDIM; dim++) {
      for (int j=0; j<N; j++) {
	Vec3 r_j = PathData.Path(slice, first+j);
	r_j[dim] += eps;
	vals.reference(Matrix(j,Range::all()));
	BandSplines(r_j[0], r_j[1], r_j[2], vals);
      }
      plus[dim] = Determinant (Matrix);

      for (int j=0; j<N; j++) {
	Vec3 r_j = PathData.Path(slice, first+j);
	r_j[dim] -= eps;
	vals.reference(Matrix(j,Range::all()));
	BandSplines(r_j[0], r_j[1], r_j[2], vals);
      }
      minus[dim] = Determinant (Matrix);
    }
    Gradient(i) = (plus-minus)/(2.0*eps);
  }
  return det;
}


void
GroundStateNodalActionClass::UpdateBands()
{
  SpeciesClass& ionSpecies = PathData.Path.Species(IonSpeciesNum);
  int first = ionSpecies.FirstPtcl;
  for (int i=0; i<NumIons; i++)
    Rions(i) = PathData.Path(0,i+first);
  System->SetIons (Rions);
  System->DiagonalizeH();
  // Now, make bands real and put into splines
  Array<double,4> data(xGrid.NumPoints, yGrid.NumPoints, zGrid.NumPoints, NumBands);
  for (int band=0; band<NumBands; band++) {
    System->SetRealSpaceBandNum(band);
    complex<double> c0 = System->RealSpaceBand(xGrid.NumPoints/2, yGrid.NumPoints/2, zGrid.NumPoints/2);
    double phi = -atan2 (c0.imag(), c0.real());
    complex<double> c(cos(phi), sin(phi));
    for (int ix=0; ix<xGrid.NumPoints-1; ix++)
      for (int iy=0; iy<yGrid.NumPoints-1; iy++)
	for (int iz=0; iz<zGrid.NumPoints-1; iz++)
	  data(ix,iy,iz,band) = real(c*System->RealSpaceBand(ix,iy,iz));
  }
  MakePeriodic (data);
  BandSplines.Init (&xGrid, &yGrid, &zGrid, data);
}


bool
GroundStateNodalActionClass::IsPositive (int slice)
{
  if (IonsHaveMoved()) 
    UpdateBands();

  double upDet, downDet;
  SpeciesClass& upSpecies = PathData.Path.Species(UpSpeciesNum);
  SpeciesClass& downSpecies = PathData.Path.Species(DownSpeciesNum);
  int first = upSpecies.FirstPtcl;
  int N = upSpecies.NumParticles;

  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = PathData.Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }
  upDet = Determinant (Matrix);

  if (upDet > 0.0) {
    first = downSpecies.FirstPtcl; 
    for (int j=0; j<N; j++) {
      Vec3 r_j = PathData.Path(slice, first+j);
      vals.reference(Matrix(j,Range::all()));
      BandSplines(r_j[0], r_j[1], r_j[2], vals);
    }
    downDet = Determinant (Matrix);
    return (downDet > 0.0);
  }
  else
    return false;
  
}
