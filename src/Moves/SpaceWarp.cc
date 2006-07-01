#include "SpaceWarp.h"

void
SpaceWarpClass::PutInBox(Vec3 &r)
{
  r[0] -= Box[0]*floor(r[0]*BoxInv[0]+0.5);
  r[1] -= Box[1]*floor(r[1]*BoxInv[1]+0.5);
  r[2] -= Box[2]*floor(r[2]*BoxInv[2]+0.5);
}

void
SpaceWarpClass::Set(const Array<Vec3,1> &rions, 
		    const Array<Vec3,1> &delta,
		    Vec3 box)
{
  Box = box;
  BoxInv = Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  N = rions.size();
  assert (delta.size() == N);
  if (Rions.size() != N) {
    Rions.resize(N);
    DeltaRions.resize(N);
    Weights.resize(N);
    Wgr.resize(N);
    rhat.resize(N);
  }
  Rions = rions;
  DeltaRions = delta;
}

Vec3
SpaceWarpClass::ForwardWarp(Vec3 r)
{
  Vec3 disp;
  double dist, distInv, totalWeight;
  for (int i=0; i<N; i++) {
    disp = r - Rions(i);
    PutInBox(disp);
    dist = sqrt(dot(disp,disp));
    distInv = 1.0/dist;
    Weights(i) = phi(distInv);
    totalWeight += Weights(i);
  }
  Weights = (1.0/totalWeight)*Weights;
  for (int i=0; i<N; i++)
    r += Weights(i)*DeltaRions(i);
  return r;
}


Vec3
SpaceWarpClass::ForwardWarp(Vec3 r,
			    TinyMatrix<double,3,3> &jmat)
{
  Vec3 disp;
  double dist, distInv, totalWeight;
  for (int i=0; i<N; i++) {
    disp = r - Rions(i);
    PutInBox(disp);
    dist = sqrt(dot(disp,disp));
    double distInv = 1.0/dist;
    Weights(i) = phi(distInv);
    totalWeight += Weights(i);
    g(i) = d_ln_phi(distInv);
    rhat(i) = distInv * disp;
  }
  Weights = (1.0/totalWeight)*Weights;
  for (int i=0; i<N; i++)
    r += Weights(i)*DeltaRions(i);
  /// Now calculate Jacobian
  jmat = 1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;
  Vec3 sumwgr (0.0, 0.0, 0.0);
  for (int j=0; j<N; j++) {
    Wgr(j) = Weights(j)*g(j)*rhat(j);
    sumwgr += Wgr(j);
  }
  for (int i=0; i<N; i++) {
    Vec3 dwdr = (Wgr(i) - Weights(i)*sumwgr);
    jmat(0,0) += DeltaRions(i)[0] * dwdr[0];
    jmat(0,1) += DeltaRions(i)[0] * dwdr[1];
    jmat(0,2) += DeltaRions(i)[0] * dwdr[2];
    jmat(1,0) += DeltaRions(i)[1] * dwdr[0];
    jmat(1,1) += DeltaRions(i)[1] * dwdr[1];
    jmat(1,2) += DeltaRions(i)[1] * dwdr[2];
    jmat(2,0) += DeltaRions(i)[2] * dwdr[0];
    jmat(2,1) += DeltaRions(i)[2] * dwdr[1];
    jmat(2,2) += DeltaRions(i)[2] * dwdr[2];
  }
  return r;
}


