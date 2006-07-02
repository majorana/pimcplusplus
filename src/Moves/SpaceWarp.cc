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
		    Vec3 box, bool doForward)
{
  Box = box;
  DoForward = doForward;
  BoxInv = Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  N = rions.size();
  assert (delta.size() == N);
  if (Rions.size() != N) {
    Rions.resize(N);
    DeltaRions.resize(N);
    Weights.resize(N);
    Wgr.resize(N);
    rhat.resize(N);
    g.resize(N);
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
  Vec3 wr;
  wr = r;
  double dist, distInv, totalWeight;
  totalWeight = 0.0;
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
    wr += Weights(i)*DeltaRions(i);
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
  return wr;
}

/// Try to find the r that would have warped me to this position. 
Vec3
SpaceWarpClass::ReverseWarp (Vec3 rp, Mat3 &jRev)
{
  Vec3 rtrial, wrtrial;
  Mat3 jForw, jForwInv;
//   // First, apply the forward warp
//   rtrial = ForwardWarp (rp, jForw);
//   // Invert the direction of warp.  This gives us our first guess for
//   // the inverse.
//   rtrial = rp - (rtrial - rp);
  
  // Initialize my guess for r with rp
  rtrial = rp;

  wrtrial = ForwardWarp (rtrial, jForw);
  do {
    jForwInv = Inverse (jForw);
    // Solve for rtrial
    rtrial = jForwInv*(rp-wrtrial) + rtrial;
    // Now calculate the forward warp from that position
    wrtrial = ForwardWarp (rtrial, jForw);
  } while (dot(rp-wrtrial,rp-wrtrial)>1.0e-30);
  jRev = Inverse(jForw);
  return rtrial;
}

void
SpaceWarpClass::SimilarTriangles 
(const Vec3 &r0,  const Vec3 &r1,  const Vec3 &r2,
 const Vec3 &r0p,       Vec3 &r1p, const Vec3 &r2p,
 double &alpha, double &beta, double &bratio)
 
{
  Vec3 a, b, bp;
  Vec3 hhat, hhatp;
  double norm, bmag, bpmag, h, hp;

  b = r2-r0; 
  bp = r2p-r0p;
  a = r1-r0;
  bmag = sqrt(dot(b,b));
  bpmag = sqrt(dot(bp,bp));
  bratio = bpmag/bmag;
  
  alpha = dot(a,b)/bmag;
  beta  = bmag - alpha;

  hhat = cross(cross(b,a),b);
  norm = 1.0/sqrt(dot(hhat,hhat));
  hhat = norm*hhat;
  h = dot (a, hhat);
  hp = bratio * h;

  hhatp = cross(cross(bp,hhat),bp);
  norm = 1.0/sqrt(dot(hhatp, hhatp));
  hhatp = norm *hhatp;

  if (dot(hhatp,hhat)< 0.0)
    hhatp = -1.0*hhatp;

  r1p = r0p + alpha/bmag *bp + hp*hhatp;

}
  
