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
  return (det/sqrt(gradMag));
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
GroundStateClass::d_dBeta (int slice1, int slice2, int level,
				      int speciesNum)
{
  return 0.0;
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
  yGrid.Init (-0.5*box[1], 0.5*box[1], nx);
  zGrid.Init (-0.5*box[2], 0.5*box[2], nx);
  Array<double,4> initData(nx,ny,nz,NumBands);
  initData = 0.0;
  BandSplines.Init (&xGrid, &yGrid, &zGrid, initData, true);
  UpdateBands();

}

double
GroundStateClass::Det(int slice, int speciesNum)
{
  SpeciesClass &species = Path.Species(speciesNum);
  int N = species.NumParticles;
  int first = species.FirstPtcl;
  assert (N == Matrix.rows());

  // First, fill up determinant matrix
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }

  // Compute determinant and cofactors
  Cofactors = Matrix;
  double det = Determinant (Cofactors);

  return det;
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
      Gradient(j) += Cofactors(i,j)*GradMat(i,j);
  }
  cerr << "Analytic gradient = " << Gradient << endl;
  GradientDetFD(slice, speciesNum);
  cerr << "FD gradient = " << Gradient << endl;
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
  Array<double,1> vals;
  for (int j=0; j<N; j++) {
    Vec3 r_j = Path(slice, first+j);
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }

  double det = Determinant (Matrix);

  for (int i=0; i<N; i++) {
    Vec3 plus, minus;
    plus = 0.0; minus = 0.0;
    for (int dim=0; dim<NDIM; dim++) {
      for (int j=0; j<N; j++) {
	Vec3 r_j = Path(slice, first+j);
	r_j[dim] += eps;
	vals.reference(Matrix(j,Range::all()));
	BandSplines(r_j[0], r_j[1], r_j[2], vals);
      }
      plus[dim] = Determinant (Matrix);

      for (int j=0; j<N; j++) {
	Vec3 r_j = Path(slice, first+j);
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
GroundStateClass::UpdateBands()
{
  SpeciesClass& ionSpecies = Path.Species(IonSpeciesNum);
  int first = ionSpecies.FirstPtcl;
  for (int i=0; i<NumIons; i++)
    Rions(i) = Path(0,i+first);
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
GroundStateClass::IsPositive (int slice, int speciesNum)
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
    vals.reference(Matrix(j,Range::all()));
    BandSplines(r_j[0], r_j[1], r_j[2], vals);
  }
  double det = Determinant (Matrix);
  cerr << "slice=" << slice 
       << " species=" << speciesNum << " det = " << det << endl;
  return (det > 0.0);  
}


bool
GroundStateNodalActionClass::IsPositive (int slice)
{
  return GroundState.IsPositive(slice, SpeciesNum);
}

double
GroundStateNodalActionClass::Action (int slice1, int slice2, 
				     const Array<int,1> &activeParticles,
				     int level)
{
  return GroundState.Action (slice1, slice2, activeParticles, level,
			     SpeciesNum);
}

double
GroundStateNodalActionClass::d_dBeta(int slice1, int slice2, int level)
{
  return GroundState.d_dBeta (slice1, slice2, level, SpeciesNum);
}
