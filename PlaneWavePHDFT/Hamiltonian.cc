#include "Hamiltonian.h"
#include "../PH/kSpacePH.h"

void GVecsClass::Set (Vec3 box, double kcut)
{
  Box = box;
  kCut = kcut;
  kBox[0]=2.0*M_PI/box[0]; kBox[1]=2.0*M_PI/box[1]; kBox[2]=2.0*M_PI/box[2];

  int maxX, maxY, maxZ;
  maxX = (int) ceil(kcut/kBox[0]);
  maxY = (int) ceil(kcut/kBox[1]);
  maxZ = (int) ceil(kcut/kBox[2]);

  // The FFT box must be twice the size as the maximum G in each direction.
  Nx = 4*maxX+1;
  Ny = 4*maxY+1;
  Nz = 4*maxZ+1;

  ////////////////////////////////////////////
  // First, set up G-vectors and difference //
  ////////////////////////////////////////////
  int numG = 0;
  Vec3 G;
  /// First, count k-vectors
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut))
	  numG++;
      }
    }
  }
  GVecs.resize(numG);
  Indices.resize(numG);
  numG = 0;
  /// Now assign
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut)) {
	  GVecs(numG) = G;
	  Indices(numG) = Int3(ix,iy,iz);
	  numG++;
	}
      }
    }
  }
  cerr << "Using " << numG << " G-vectors.\n";
  ////////////////////////////////////////////
  // Now, set up G-vector differences.      //
  ////////////////////////////////////////////
  numG = 0;
  /// First, count k-vectors
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (4.0*kcut*kcut))
	  numG++;
      }
    }
  }
  GDiff.resize(numG);
  IDiff.resize(numG);
  numG = 0;
  /// Now assign
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (4.0*kcut*kcut)) {
	  GDiff(numG) = G;
	  IDiff(numG) = Int3(ix,iy,iz);
	  numG++;
	}
      }
    }
  }
}


void GVecsClass::GetFFTBoxSize(int &nx, int &ny, int &nz)
{
  nx=Nx; ny=Ny; nz=Nz;
}

void
FFTBox::PutkVec (const zVec &c)
{
  kBox = 0.0;
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.Index(n)[0]+Nx)%Nx;
      int j = (GVecs.Index(n)[1]+Ny)%Ny;
      int k = (GVecs.Index(n)[2]+Nz)%Nz;
      kBox(i,j,k) = c(n);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
      int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
      int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
      kBox(i,j,k) = c(n);
    }
  else {
    cerr << "Incommensurate dimensions in PutkVec.\n";
    abort();
  }
}

void 
FFTBox::GetkVec (zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.Index(n)[0]+Nx)%Nx;
      int j = (GVecs.Index(n)[1]+Ny)%Ny;
      int k = (GVecs.Index(n)[2]+Nz)%Nz;
      c(n) = kBox(i,j,k);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
      int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
      int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
      c(n) = kBox(i,j,k);
    }
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


void
FFTBox::AddFromVec (const zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.Index(n)[0]+Nx)%Nx;
      int j = (GVecs.Index(n)[1]+Ny)%Ny;
      int k = (GVecs.Index(n)[2]+Nz)%Nz;
      kBox(i,j,k) += c(n);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
      int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
      int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
      kBox(i,j,k) += c(n);
    }
  else {
    cerr << "Incommensurate dimensions in AddFromVec.\n";
    abort();
  }
}

void 
FFTBox::AddToVec (zVec &c)
{
  if (c.size() == GVecs.size()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.Index(n)[0]+Nx)%Nx;
      int j = (GVecs.Index(n)[1]+Ny)%Ny;
      int k = (GVecs.Index(n)[2]+Nz)%Nz;
      c(n) += kBox(i,j,k);
    }
  else if (c.size() == GVecs.DeltaSize()) 
    for (int n=0; n<c.size(); n++) {
      int i = (GVecs.DeltaI(n)[0]+Nx)%Nx;
      int j = (GVecs.DeltaI(n)[1]+Ny)%Ny;
      int k = (GVecs.DeltaI(n)[2]+Nz)%Nz;
      c(n) += kBox(i,j,k);
    }
  else {
    cerr << "Incommensurate dimensions in GetkVec.\n";
    abort();
  }
}


void KineticClass::Setup()
{
  halfG2.resize(GVecs.size());
  for (int i=0; i<GVecs.size(); i++)
    halfG2(i) = 0.5*dot(GVecs(i),GVecs(i));
  IsSetup = true;
}

void KineticClass::Apply(const zVec &c, zVec &Kc)
{
  if (!IsSetup)
    Setup();
  for (int i=0; i<c.size(); i++)
    Kc(i) += halfG2(i)*c(i);
}


void CoulombClass::Apply(const zVec &c, zVec &Vc)
{
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<c.size(); i++) {
    for (int j=0; j<i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      Vc(i) += -4.0*volInv*M_PI*Z/dot(diff,diff)*c(j);
    }
    for (int j=i+1; j<c.size(); j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      Vc(i) += -4.0*volInv*M_PI*Z/dot(diff,diff)*c(j);
    }
  }
}


void CoulombFFTClass::Setup()
{
  FFT.Setup();
  Vec3 dG;
  Int3 dI;
  //double prefact = -4.0*M_PI*Z;// /GVecs.GetBoxVol();
  double prefact = -4.0*M_PI*Z/GVecs.GetBoxVol();
  FFT.kBox = 0.0;
//   for (int m=0; m<GVecs.size(); m++)
//     for (int n=0; n<GVecs.size(); n++) {
//       dG = GVecs(m)-GVecs(n);
//       dI = GVecs.Index(n) - GVecs.Index(m);
//       double val = (m == n) ? 0.0 : prefact/dot(dG,dG);
//       int i = (dI[0] + FFT.Nx) % FFT.Nx;
//       int j = (dI[1] + FFT.Ny) % FFT.Ny;
//       int k = (dI[2] + FFT.Nz) % FFT.Nz;
//       FFT.kBox(i,j,k) = val;
//     }
  Vec3 Gdiff;
  int Nx=FFT.Nx; int Ny=FFT.Ny; int Nz=FFT.Nz;
  cerr << "Nx=" << Nx << " Ny=" << Ny << " Nz=" << Nz << endl;
  int mx = (Nx-1)/2;
  int my = (Ny-1)/2;
  int mz = (Nz-1)/2;
  double cut = 4.0*GVecs.GetkCut() * GVecs.GetkCut();
  for (int ix=-mx; ix<=mx; ix++) {
    int nx = (ix+Nx)%Nx;
    Gdiff[0] = GVecs.GetkBox()[0]*ix;
    for (int iy=-my; iy<=my; iy++) {
      int ny = (iy+Ny)%Ny;
      Gdiff[1] = GVecs.GetkBox()[1]*iy;
      for (int iz=-mz; iz<=mz; iz++) {
	int nz = (iz+Nz)%Nz;
	Gdiff[2] = GVecs.GetkBox()[2]*iz;
	double prod = dot(Gdiff,Gdiff);
	if (prod < cut) {
	  //cerr << "oldval = " << FFT.kBox(nx,ny,nz) <<endl;
	  FFT.kBox(nx,ny,nz) = prefact/prod;
	  //cerr << "newval = " << FFT.kBox(nx,ny,nz) <<endl;
	}
      }
    }
  }
  FFT.kBox(0,0,0) = 0.0;
    
  FFT.k2r();
  FILE *fout = fopen ("Vr.dat", "w");
  for (int i=0; i<FFT.Nx; i++) {
    for (int j=0; j<FFT.Ny; j++)
      fprintf (fout, "%1.15e ", FFT.rBox(i,j,0).real());
    fprintf (fout, "\n");
  }
  fclose (fout);
  Vr.resize(FFT.Nx, FFT.Ny, FFT.Nz);
  Vr = FFT.rBox;
  IsSetup = true;
  cerr << "Finished CoulombFFT setup.\n";
}

void 
CoulombFFTClass::Apply(const zVec &c, zVec &Vc)
{
  if (!IsSetup)
    Setup();

  double Ninv = 1.0/(FFT.Nx*FFT.Ny*FFT.Nz);
  // First, put c into FFTbox
  FFT.PutkVec (c);
  // Now, transform to real space
  FFT.k2r();
  FFT.rBox *= sqrt(Ninv);
  // Now, multiply by V
  FFT.rBox *= Vr;
//   for (int i=0; i<FFT.Nx; i++)
//     for (int j=0; j<FFT.Ny; j++)
//       for (int k=0; k<FFT.Nz; k++)
// 	FFT.rBox(i,j,k) *= Vr(i,j,k); 

  // Transform back
  FFT.r2k();
  FFT.kBox *= sqrt(Ninv);
  // And put into Vc
  FFT.AddToVec (Vc);
}


void PHPotClass::Setup()
{
  kPH.CalcTailCoefs (30.0, 80.0);
  VGGp.resize (GVecs.size(), GVecs.size());
  Vec3 k (0.0, 0.0, 0.0);
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<GVecs.size(); i++)
    for (int j=i; j<GVecs.size(); j++) {
      VGGp (i,j) = kPH.V(k, GVecs(i), GVecs(j))*volInv;
      VGGp (j,i) = conj (VGGp(i,j));
    }

  IsSetup = true;
}

void PHPotClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
      
  for (int i=0; i<c.size(); i++)
    for (int j=0; j<c.size(); j++)
      Hc(i) += VGGp (i,j) * c(j);
}


void Hamiltonian::Apply(const zVec &c, zVec &Hc)
{
  Hc = 0.0;
  Kinetic.Apply (c, Hc);
  //Coulomb.Apply (c, Hc);
  //CoulombFFT.Apply (c, Hc);
  PH.Apply (c, Hc);
}


