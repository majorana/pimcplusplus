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
GroundStateNodalActionClass::Action (int slice1, int slice2, 
				     const Array<int,1> &activeParticles, 
				     int level)
{
  double lambdaTauInv = 1.0/(Path.tau*Path.Species(UpSpeciesNum).lambda);

  // The nodal action should only be used at level = 0;
  assert (level == 0);
  bool doUp   = false;
  bool doDown = false;
  if (IonsHaveMoved()) {
    
  }
  for (int i=0; i<activeParticles.size(); i++) 
    if (Path.ParticleSpeciesNum(activeParticles(i)) == UpSpeciesNum)
      doUp = true;
    else if (Path.ParticleSpeciesNum(activeParticles(i)) == DownSpeciesNum)
      doDown = true;

  double action = 0.0;
  if (doUp) {
    for (int slice=slice1; slice <= slice2; slice++) {
      double det = GradientDet (slice, UpSpeciesNum);
      double gradMag = 0.0;
      for (int i=0; i<Path.Species(UpSpeciesNum).NumParticles; i++)
	gradMag += dot (Gradient(i), Gradient(i));
      // Newton Raphson approximation to distance to node
      UpDists(slice) = det/sqrt(gradMag);
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
  
  int numUp   =  PathData.Path.Species(UpSpeciesNum).NumParticles;
  int numDown = (DownSpeciesNum != -1) ? 
    PathData.Path.Species(DownSpeciesNum).NumParticles : -1;

  if (DownSpeciesNum != -1)
    assert (numUp == numDown);
  
  Workspace.resize (DetCofactorsWorksize(numUp));
  Matrix.resize(numUp, numUp);
  Cofactors.resize(numUp, numUp);
  Gradient.resize(numUp);
  TempGrad.resize(numUp);
  UpDists.resize(PathData.NumTimeSlices());
  DownDists.resize(PathData.NumTimeSlices());

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
    BandSplines.Grad(r_i[0], r_i[1], r_i[2], TempGrad);
    for (int j=0; j<N; j++)
      Gradient(j) += Cofactors(i,j)*TempGrad(j);
  }
  return det;
}
