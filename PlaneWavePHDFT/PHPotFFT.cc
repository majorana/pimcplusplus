#include "PHPotFFT.h"


void
PHPotFFTClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  if (!IsSetup)
    Setup();
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<vmat.rows(); i++) 
    for (int j=0; j<=i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      vmat (i,j) = s*kPH.V(kPoint, GVecs(i), GVecs(j))*volInv;
      vmat (j,i) = conj(vmat(i,j));
    }
//   for (int i=0; i<vmat.rows(); i++)
//     cerr << "vmat(" << i << "," << i << ") = " << vmat(i,i) <<
//     endl;
//  cerr << "vmat = " << vmat << endl;
}


void 
PHPotFFTClass::SetupkPotentials()
{
  kPH.CalcTailCoefs (30.0, 60.0);

  double volInv = 1.0/GVecs.GetBoxVol();
  // Setup V and F tensors in k-space
  //  double gMag, lastMag2;
  double gMag, lastMag2;
  lastMag2 = -1.0;
  double a, bPerp, bPar, V;
  int numCalls = 0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    double gMag2 = dot(GVecs.DeltaG(i), GVecs.DeltaG(i));
    if (fabs(gMag2-lastMag2) > 1.0e-12) {
      lastMag2 = gMag2;
      gMag  = sqrt(lastMag2);
      kPH.GetVals(gMag, a, bPerp, bPar, V);
      numCalls++;
    }
    Vk(i) = volInv * V;
    Fk(i) = volInv * kPH.Ftensor (GVecs.DeltaG(i), a, bPerp, bPar);
  }
}

void 
PHPotFFTClass::SetuprPotentials()
{
  cFFT.kBox = complex<FFT_FLOAT>(0.0, 0.0);
  complex<FFT_FLOAT> z(0.0, 0.0);
  TinyMatrix<complex<FFT_FLOAT>,3,3> Mat0;
  Mat0(0,0)=z; Mat0(0,1)=z; Mat0(0,2)=z; 
  Mat0(1,0)=z; Mat0(1,1)=z; Mat0(1,2)=z;
  Mat0(2,0)=z; Mat0(2,1)=z; Mat0(2,2)=z;
  MatFFT.kBox = Mat0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Int3 I = GVecs.DeltaI(i);
    cFFT.kBox(I)   = StructureFactor(i)*Vk(i);
#ifdef FFT_USE_SINGLE
    MatFFT.kBox(I) = conv(StructureFactor(i)*Fk(i));
#else
    MatFFT.kBox(I) = StructureFactor(i)*Fk(i);
#endif
  }
  cFFT.k2r();
  Vr = cFFT.rBox;
  MatFFT.k2r();
  Fr = MatFFT.rBox;
}

void
PHPotFFTClass::SetIons(const Array<Vec3,1> &rions)
{
  // Calculate the structure factor
  VionBase::SetIons(rions);
  if (IsSetup)
    SetuprPotentials();
}


void 
PHPotFFTClass::Setup()
{
  cerr << "Setting up FFT boxes:\n";
  VecFFT.Setup();
  MatFFT.Setup();
  cerr << "done.\n";

  kPH.CalcTailCoefs (15.0, 60.0);
  int nx, ny, nz;
  MatFFT.GetDims(nx,ny,nz);
  Fr.resize(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Gc.resize(GVecs.size());
  Vc.resize(GVecs.size());
  Vk.resize(GVecs.DeltaSize());
  Fk.resize(GVecs.DeltaSize());
//   StructureFactor.resize(GVecs.DeltaSize());
//   StructureFactor = complex<double>(1.0, 0.0);

  SetupkPotentials();
  SetuprPotentials();

  IsSetup = true;
}
  
void
PHPotFFTClass::Setk (Vec3 k)
{
  kPoint = k;
  VecFFT.Setup();
  MatFFT.Setup();
  
  int nx, ny, nz;
  MatFFT.GetDims(nx,ny,nz);
  Fr.resize(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Gc.resize(GVecs.size());
  Vc.resize(GVecs.size());
  Vk.resize(GVecs.DeltaSize());
  Fk.resize(GVecs.DeltaSize());

  SetupkPotentials();
  SetIons(Rions);
}


void 
PHPotFFTClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);
  double nInv = 1.0/(double)(nx*ny*nz);

  /////////////////////////
  // Pseudo-kinetic part //
  /////////////////////////
  // First, build G*c vector in VecFFT  
  for (int i=0; i<GVecs.size(); i++) 
    Gc(i) = c(i)*(kPoint+GVecs(i));
  VecFFT.PutkVec (Gc);
  VecFFT.k2r();
  
  // Now, muliply by F(r)
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++) {
	register complex<FFT_FLOAT> vx = VecFFT.rBox(ix,iy,iz)[0];
	register complex<FFT_FLOAT> vy = VecFFT.rBox(ix,iy,iz)[1];
	register complex<FFT_FLOAT> vz = VecFFT.rBox(ix,iy,iz)[2];
	register complex<FFT_FLOAT> Fvx, Fvy, Fvz;
	Fvx  = Fr(ix,iy,iz)(0,0) * vx;
	Fvy  = Fr(ix,iy,iz)(1,0) * vx;
	Fvz  = Fr(ix,iy,iz)(2,0) * vx;
	Fvx += Fr(ix,iy,iz)(0,1) * vy;
	Fvy += Fr(ix,iy,iz)(1,1) * vy;
	Fvz += Fr(ix,iy,iz)(2,1) * vy;
	Fvx += Fr(ix,iy,iz)(0,2) * vz;
	Fvy += Fr(ix,iy,iz)(1,2) * vz;
	Fvz += Fr(ix,iy,iz)(2,2) * vz;
	VecFFT.rBox(ix,iy,iz) =
	TinyVector<complex<FFT_FLOAT>,3>(Fvx,Fvy,Fvz);

//	VecFFT.rBox(ix,iy,iz) = Fr(ix,iy,iz)*VecFFT.rBox(ix,iy,iz);

// 	VecFFT.rBox(ix,iy,iz)[0] = 
// 	  Fr(ix,iy,iz)(0,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(0,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(0,2)*VecFFT.rBox(ix,iy,iz)[2];
// 	VecFFT.rBox(ix,iy,iz)[1] = 
// 	  Fr(ix,iy,iz)(1,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(1,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(1,2)*VecFFT.rBox(ix,iy,iz)[2];
// 	VecFFT.rBox(ix,iy,iz)[2] = 
// 	  Fr(ix,iy,iz)(2,0)*VecFFT.rBox(ix,iy,iz)[0] + 
// 	  Fr(ix,iy,iz)(2,1)*VecFFT.rBox(ix,iy,iz)[1] + 
// 	  Fr(ix,iy,iz)(2,2)*VecFFT.rBox(ix,iy,iz)[2];
      }
  // Transform back;
  VecFFT.r2k();
  // Rename for clarity
  zVecVec& FGc = Gc;
  VecFFT.GetkVec (FGc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += 0.5*nInv*dot(GVecs(i)+kPoint, FGc(i));

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();
  // Multiply by V
  cFFT.rBox *= Vr;
  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += nInv*Vc(i);
}

