#include "FFTBox.h"

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


////////////////////////////////
// FFTVecBox Member Functions //
////////////////////////////////
void
FFTVecBox::PutkVec (const zVecVec &c)
{
  complex<double> zero(0.0, 0.0);
  cVec3 zero3(zero, zero, zero);
  kBox = zero3;
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
FFTVecBox::GetkVec (zVecVec &c)
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
FFTVecBox::AddFromVec (const zVecVec &c)
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
FFTVecBox::AddToVec (zVecVec &c)
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


////////////////////////////////
// FFTMatBox Member Functions //
////////////////////////////////
void
FFTMatBox::PutkVec (const zMatVec &c)
{
  cMat3 zero;
  zero(0,0)=0.0; zero(0,1)=0.0; zero(0,2)= 0.0;
  zero(1,0)=0.0; zero(1,1)=0.0; zero(1,2)= 0.0;
  zero(2,0)=0.0; zero(2,1)=0.0; zero(2,2)= 0.0;
  kBox = zero;
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
FFTMatBox::GetkVec (zMatVec &c)
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
FFTMatBox::AddFromVec (const zMatVec &c)
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
FFTMatBox::AddToVec (zMatVec &c)
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
