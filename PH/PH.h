#ifndef PH_H
#define PH_H

#include "../Integration/Integrate.h"
#include "../Splines/Grid.h"
#include "Chebyshev.h"
#include "../Splines/CubicSpline.h"
#include "../IO/InputOutput.h"

enum PHType {PH_NONE, PH_CHEBYSHEV, PH_CUBIC, PH_CUBICXC, PH_NUCLEAR};

class FullCorePotential
{
private:
  int GridInitialized;
public:
  char FileName[1000];

  CubicSpline V;
  scalar Z;

  inline scalar operator()(scalar r)
  {
    if (r < V.grid->End)
      return (-Z/r + V(r));
    else
      return (-Z/r);
  }
  inline scalar Deriv(scalar r)
  {
    if (r < V.grid->End)
      return (Z/(r*r) + V.Deriv(r));
    else      
      return (Z/(r*r));
  }
  void Read(char *FName);
  void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("Z", Z);
    outSection.OpenSection("V");
    V.Write(outSection);
    outSection.CloseSection();
  }
  void Read  (InputSectionClass &inSection)
  {
    assert(inSection.ReadVar("Z", Z));
    assert(inSection.OpenSection("V"));
    V.Read (inSection);
    inSection.CloseSection();
  }

  FullCorePotential()
  {
    GridInitialized = 0;
  }
  ~FullCorePotential()
  {
    if (GridInitialized)
      delete V.grid;
  }
};




class PseudoHamiltonian
{
public:
  scalar CoreRadius;
  scalar Z, Zion;  // Zion is the unscreened pseudohamiltonian charge
  double Amin, Amax, Bmin, Bmax, Vmin, Vmax;

  FullCorePotential *FullCoreV;
  CubicSpline VHXC;  // The Hartree, exchange, and correlation potentials
  scalar NumElecs;
  int UseVHXC;

  inline void UpdateVHXC (CubicSpline NewVHXC, scalar Nelecs)
  {
    NumElecs = Nelecs;
    VHXC = NewVHXC;
    UseVHXC = 1;
  }

  virtual PHType Type()
  {
    cerr << "PH base class has not type.\n";
    exit(1);
    return (PH_NONE);
  }

  virtual int NumParams()
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0);
  }

  virtual scalar &Params(int i)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    static scalar r = 0.0;
    return(r);
  }
  
  virtual scalar Params(int i) const
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		   scalar &dAdr)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }

  virtual scalar d2Adr2(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    return (0.0);
  }

  virtual scalar dVdr(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual scalar    V(scalar r)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
    exit(1);
    return(0.0);
  }
  virtual void Write(char *FileName)
  {
    FILE *fout;
    if ((fout = fopen (FileName, "w")) == NULL)
      {
	cerr << "Can't open " << FileName << " for writing.  Exitting.\n";
	exit(1);
      }
    Write (fout);
  }
  virtual void Write(FILE *fout)

  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }
  virtual void Read (FILE *fin)
  {
    cerr << "Should never get to PseudoHamiltonian base class.\n";
  }

  virtual void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("CoreRadius", CoreRadius);
    outSection.WriteVar ("Z", Z);
    outSection.WriteVar ("Zion", Zion);
    outSection.WriteVar ("NumElecs", NumElecs);
    outSection.WriteVar ("Amin", Amin);
    outSection.WriteVar ("Amax", Amax);
    outSection.WriteVar ("Bmin", Bmin);
    outSection.WriteVar ("Bmax", Bmax);
    outSection.WriteVar ("Vmin", Vmin);
    outSection.WriteVar ("Vmax", Vmax);
    outSection.WriteVar ("UseVHXC", UseVHXC);
    if (UseVHXC) {
      outSection.OpenSection("VHXC");
      VHXC.Write (outSection);
      outSection.CloseSection();
    }
    outSection.OpenSection ("FullCoreV");
    FullCoreV->Write(outSection);
    outSection.CloseSection();
  }

  virtual bool Read (InputSectionClass &inSection)
  {
    assert (inSection.ReadVar ("CoreRadius", CoreRadius));
    assert (inSection.ReadVar ("Z", Z));
    assert (inSection.ReadVar ("Zion", Zion));
    assert (inSection.ReadVar ("NumElecs", NumElecs));
    assert (inSection.ReadVar ("Amin", Amin));
    assert (inSection.ReadVar ("Amax", Amax));
    assert (inSection.ReadVar ("Bmin", Bmin));
    assert (inSection.ReadVar ("Bmax", Bmax));
    assert (inSection.ReadVar ("Vmin", Vmin));
    assert (inSection.ReadVar ("Vmax", Vmax));
    assert (inSection.ReadVar ("UseVHXC", UseVHXC));
    if (UseVHXC){
      assert(inSection.OpenSection("VHXC"));
      VHXC.Read (inSection);
      inSection.CloseSection();
    }
    assert(inSection.OpenSection ("FullCoreV"));
    FullCoreV = new FullCorePotential;
    FullCoreV->Read (inSection);
    inSection.CloseSection();
  } 



  PseudoHamiltonian()
  {
    UseVHXC = 0;
  }
};

///////////////////////////////////////////////////////////////////
//                       Cubic Spline PH                         //
///////////////////////////////////////////////////////////////////

class PH_CubicSpline : public PseudoHamiltonian
{
private:
  int UpDateB;
public:
  CubicSpline PA, PB, Vfunc;
  LinearGrid Agrid, Bgrid, Vgrid;
  
  PHType Type()
  {
    return PH_CUBIC;
  }

  int NumParams()
  {
    return (PA.NumParams-1 + PB.NumParams-2 + Vfunc.NumParams-1);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      UpDateB = 1;
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

  scalar Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	scalar pa = PA(r);
	scalar pb = PB(r);
	A = Amin + pa * pa;
	B = Bmin + pb * pb;
	V = Vfunc(r);
	dAdr = 2.0*pa * PA.Deriv(r);
      }
    else
      {
	A = 1.0;
	B = 1.0;
	V = (*FullCoreV)(r);
	dAdr = 0.0;
      }
  }

  scalar V(scalar r)
  {
    if (r<CoreRadius)
      return (Vfunc(r));
    else
      return ((*FullCoreV)(r));
  }

  scalar dVdr(scalar r)
  {
    if (r<CoreRadius)
      return (Vfunc.Deriv(r));
    else
      return (FullCoreV->Deriv(r));
  }

  scalar d2Adr2(scalar r)
  {
    if (r<CoreRadius)
      {
	scalar p = PA(r);
	scalar pp = PA.Deriv(r);
	scalar ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }

  void SetCoreRadius (scalar radius)
    {
      CoreRadius = radius;
      Agrid.Init(0.0, CoreRadius, PA.NumParams);
      Bgrid.Init(0.0, CoreRadius, PB.NumParams);
      int NumV = Vfunc.NumParams;
      Vgrid.Init(0.0, CoreRadius, NumV);
      Array<scalar,1> Vparams(NumV);
      for (int i=0; i<NumV; i++)
	Vparams(i) = Vfunc.Params(i);
      Vparams(NumV-1) = (*FullCoreV)(CoreRadius);
      Vfunc.Init (&Vgrid, Vparams, 5.0e30, FullCoreV->Deriv(CoreRadius));
    }
      

  inline void Init (scalar amin, scalar bmin, scalar vmin,
		    scalar amax, scalar bmax, scalar vmax,
		    Array<scalar,1> Aparams,
		    Array<scalar,1> Bparams,
		    Array<scalar,1> Vparams,
		    FullCorePotential *FullCorePot,
		    scalar coreradius)
  {
    CoreRadius = coreradius;
    Amin = amin; Amax = amax; 
    Bmin = bmin; Bmax = bmax;
    Vmin = vmin; Vmax = vmax;
    FullCoreV = FullCorePot;

    int NumA = Aparams.rows();
    Agrid.Init(0.0, CoreRadius, NumA+1);
    Array<scalar,1> Ainit(NumA+1);
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-amin);
    PA.Init(&Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Bgrid.Init(0.0, CoreRadius, NumB+2);
    Array<scalar,1> Binit(NumB+2);
    Binit(Range(1,NumB)) = Bparams;
    Binit(0) = sqrt(Aparams(0)*Aparams(0) + Amin - Bmin);    
    Binit(NumB+1) = sqrt(1.0 - bmin);
    PB.Init(&Bgrid, Binit, 0.0, 0.0);

    int NumV = Vparams.rows();
    Vgrid.Init(0.0, CoreRadius, NumV+1);
    Array<scalar,1> Vinit(NumV+1);
    Vinit(Range(0,NumV-1)) = Vparams;
    Vinit(NumV) = (*FullCoreV)(CoreRadius);
    Vfunc.Init(&Vgrid, Vinit, 5.0e30, FullCoreV->Deriv(CoreRadius));
  }
    
  PH_CubicSpline()
  {
    // Do nothing for now
  }

  void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "PH_CubicSpline");
    // Write common variables
    PseudoHamiltonian::Write(outSection);
    outSection.OpenSection ("Agrid");
    Agrid.Write(outSection);
    outSection.CloseSection();
    outSection.OpenSection ("Bgrid");
    Bgrid.Write(outSection);
    outSection.CloseSection();
    outSection.OpenSection ("Vgrid");
    Vgrid.Write(outSection);
    outSection.CloseSection();
    Array<double,1> Params;
    Params.resize(Agrid.NumPoints);
    for (int i=0; i<Agrid.NumPoints; i++)
      Params(i) = PA.Params(i);
    outSection.WriteVar ("PAparams", Params);

    Params.resize(Bgrid.NumPoints);
    for (int i=0; i<Bgrid.NumPoints; i++)
      Params(i) = PB.Params(i);
    outSection.WriteVar ("PBparams", Params);

    Params.resize(Vgrid.NumPoints);
    for (int i=0; i<Vgrid.NumPoints; i++)
      Params(i) = Vfunc.Params(i);
    outSection.WriteVar ("Vfuncparams", Params);
  }

  bool Read (InputSectionClass &inSection)
  {
    PseudoHamiltonian::Read(inSection);
    assert(inSection.OpenSection("Agrid"));
    Agrid.Read(inSection);
    inSection.CloseSection();
    
    assert(inSection.OpenSection("Bgrid"));
    Bgrid.Read(inSection);
    inSection.CloseSection();

    assert(inSection.OpenSection("Vgrid"));
    Vgrid.Read(inSection);
    inSection.CloseSection();
   
    Array<double,1> Params;
    assert(inSection.ReadVar ("PAparams", Params));
    PA.Init (&Agrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("PBparams", Params));
    PB.Init (&Bgrid, Params, 0.0, 0.0);
    assert(inSection.ReadVar ("Vfuncparams", Params));
    Vfunc.Init(&Vgrid, Params, 5.0e30, FullCoreV->Deriv(CoreRadius));
    return (true);
   }


  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Cubic Spline PH with Exchange-Correlation subtracted from V   //
// This is used to store a PH which has had the valance exchange //
// correlation subtracted off.  There is a dense grid for V,     //
// while the A anb B grids remain the same.                      //
///////////////////////////////////////////////////////////////////
class PH_CubicSplineXC : public PseudoHamiltonian
{
private:
  int UpDateB;
public:
  CubicSpline PA, PB, Vfunc;
  Grid *Agrid, *Bgrid, *Vgrid;

  // Nelecs is the number of electrons in the pseudo-atom;
  scalar NumElecs;
  scalar Amin, Bmin;

  PHType Type()
  {
    return PH_CUBICXC;
  }
  
  int NumParams()
  {
    return (PA.NumParams-1 + PB.NumParams-2 + Vfunc.NumParams-1);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      UpDateB = 1;
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

  scalar Params(int i) const
  {
    if (i < (PA.NumParams-1))
      return (PA.Params(i));
    else if (i < (PA.NumParams-1 + PB.NumParams-2))
      return (PB.Params(i-PA.NumParams + 2));
    else 
      return (Vfunc.Params(i-PA.NumParams - PB.NumParams + 3));
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    if (r<CoreRadius)
      {
	if (UpDateB)
	  {
	    PB.Params(0) = sqrt(Amin - Bmin + PA.Params(0)*PA.Params(0));
	    UpDateB = 0;
	  }
	scalar pa = PA(r);
	scalar pb = PB(r);
	A = Amin + pa * pa;
	B = Bmin + pb * pb;

	dAdr = 2.0*pa * PA.Deriv(r);
      }
    else
      {
	A = 1.0;
	B = 1.0;
	dAdr = 0.0;
      }
    if (r < Vgrid->End)
      {
	V = Vfunc(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    V += VHXC(r);
      }
    else
      V = -(Zion - NumElecs) / r;

  }

  scalar V(scalar r)
  {
    scalar v;
    if (r < Vgrid->End)
      {
	v = Vfunc(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    v += VHXC(r);
      }
    else
      v = -(Zion-NumElecs) / r;
    return (v);
  }

  scalar dVdr(scalar r)
  {
    scalar dV;
    if (r < Vgrid->End)
      {
	dV = Vfunc.Deriv(r);
	if (UseVHXC)
	  if (r < VHXC.grid->End)
	    dV += VHXC.Deriv(r);
      }
    else
      dV = (Zion-NumElecs) / (r*r);

    return (dV);
  }

  
   scalar d2Adr2(scalar r)
  {
    if (r<CoreRadius)
      {
	scalar p = PA(r);
	scalar pp = PA.Deriv(r);
	scalar ppp = PA.Deriv2(r);
	return (2.0*(p*ppp+pp*pp));
      }
    else
      return (0.0);
  }   

  inline void Init (Array<scalar,1> Aparams,
		    Array<scalar,1> Bparams,
		    Array<scalar,1> Vparams,
		    Grid *agrid, Grid *bgrid,
		    Grid *vgrid, scalar zion,
		    scalar coreradius)
  {
    UseVHXC = 0;
    NumElecs = 0;
    CoreRadius = coreradius;
    Agrid = agrid;
    Bgrid = bgrid;
    Vgrid = vgrid;
    Zion = zion;
    Amin = 0.0; 
    Bmin = 0.0;

    int NumA = Aparams.rows();
    Array<scalar,1> Ainit(NumA+1);
    if (Ainit.rows() != Agrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	exit(1);
      }
    Ainit(Range(0,NumA-1)) = Aparams;
    Ainit(NumA) = sqrt(1.0-Amin);
    PA.Init(Agrid, Ainit, 0.0, 0.0);

    int NumB = Bparams.rows();
    Array<scalar,1> Binit(NumB+2);
    if (Binit.rows() != Bgrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	exit(1);
      }
    Binit(Range(1,NumB)) = Bparams;
    Binit(0) = sqrt(Aparams(0)*Aparams(0) + Amin - Bmin);    
    Binit(NumB+1) = sqrt(1.0 - Bmin);
    PB.Init(Bgrid, Binit, 0.0, 0.0);

    int NumV = Vparams.rows();
    Array<scalar,1> Vinit(NumV+1);
    if (Vinit.rows() != Vgrid->NumPoints)
      {
	cerr << "Size mismatch in PH_CubicSplineXC.\n";
	cerr << "Vinit.rows() = " << Vinit.rows() << "\n";
	cerr << "Vgrid.NumPoints = " << Vgrid->NumPoints << "\n";
	exit(1);
      }
    Vinit(Range(0,NumV-1)) = Vparams;
    scalar rmax = Vgrid->End; 

    // Force the last point to agree with asymptotic form
    Vinit(NumV) = -Zion/rmax;
    // Make the slope match assymptotic form as well
    Vfunc.Init(Vgrid, Vinit, 5.0e30, Zion/(rmax*rmax));
  }
    
  PH_CubicSplineXC()
  {
    // Do nothing for now
  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Nuclear PH:  has no parameters.  The potential is dependent   //
// only on the nuclear charge Z, and possibly on the exchange-   //
// correlation potential.  This is used for full-core atom       //
// calculations.                                                 //
///////////////////////////////////////////////////////////////////
class PH_Nuclear : public PseudoHamiltonian
{
public:
  scalar dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = -Zion/r;
    if (UseVHXC)
      if (r < VHXC.grid->End)
	V += VHXC(r);
      else
	V = -(Zion-NumElecs)/r;
  }

  scalar V(scalar r)
  {
    scalar v = -Zion/r;
    if (UseVHXC)
      if (r < VHXC.grid->End)
	v += VHXC(r);
      else
	v = -(Zion-NumElecs) / r;
    return (v);
  }

  scalar dVdr(scalar r)
  {
    scalar dV;
    dV = Zion / (r*r);
    if (UseVHXC)
      if (r < VHXC.grid->End)
	dV += VHXC.Deriv(r);
      else
	dV = (Zion-NumElecs) / (r*r);

    return (dV);
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  bool Read (InputSectionClass &inSection)
  {
    assert(inSection.ReadVar ("Zion", Zion));
    assert(inSection.ReadVar ("UseVHXC", UseVHXC));
    if (UseVHXC)
      {
	assert(inSection.ReadVar ("NumElecs", NumElecs));
	assert(inSection.OpenSection("VHXC"));
	VHXC.Read(inSection);
	inSection.CloseSection();
      }
    return (true);
  }

  void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Nuclear");
    outSection.WriteVar ("Zion", Zion);
    outSection.WriteVar ("UseVHXC", UseVHXC);
    if (UseVHXC) {
      outSection.WriteVar ("NumElecs", NumElecs);
      outSection.OpenSection("VHXC");
      VHXC.Write(outSection);
      outSection.CloseSection();
    } 
  }

  PH_Nuclear(scalar z)
  {

    Zion = z;
  }

  PH_Nuclear()
  {
    // Do nothing for now
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



///////////////////////////////////////////////////////////////////
// FullCore PH:  has no parameters.  It is use solely to         //
// calculate the full core radial WFs for comparison with the    //
// pseudo radial WFs.                                            //
///////////////////////////////////////////////////////////////////
class PH_FullCore : public PseudoHamiltonian
{
public:
  // Zion is the unscreened pseudohamiltonian charge and Nelecs
  // is the number of electrons in the pseudo-atom;
  scalar dummy;  // To satisfy compiler warning.


  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = (*FullCoreV)(r);
  }

  scalar V(scalar r)
  {
    return ((*FullCoreV)(r));
  }

  scalar dVdr(scalar r)
  {
    return (FullCoreV->Deriv(r));
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }
  
  
  PH_FullCore(FullCorePotential *FullCorePot)
  {
    Z = FullCorePot->Z;
    FullCoreV = FullCorePot;
  }

  PH_FullCore()
  {

  }

  void Write(FILE *fout);
  void Write(char *FileName);
  void Read (FILE *fin);
  
};




///////////////////////////////////////////////////////////////////
// Zero PH:  has no parameters.  Has unit masses and V(r) = 0.   //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Zero : public PseudoHamiltonian
{
public:
  scalar dummy;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (dummy);
  }

  scalar Params(int i) const
  {
    cerr << "Nuclear PH has no parameters.\n";
    return (0.0);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = 0.0;
  }

  scalar V(scalar r)
  {
    return (0.0);
  }

  scalar dVdr(scalar r)
  {
    return (0.0);
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Zero");
  }

  PH_Zero()
  {
    // Do nothing for now
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



///////////////////////////////////////////////////////////////////
// Gaussian PH:  has two parameters.  Has unit masses and        //
// V(r) = Amp * exp(-r^2/(2.0*sigma^2)).                         //
// It is used for testing in matrix squaring.                    //
///////////////////////////////////////////////////////////////////
class PH_Gaussian : public PseudoHamiltonian
{
public:
  scalar Amp, sigma;  // To satisfy compiler warnings about no return value

  PHType Type()
  {
    return PH_NUCLEAR;
  }
  
  int NumParams()
  {
    return (0);
  }
  
  scalar &Params(int i)
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

  scalar Params(int i) const
  {
    if (i==0)
      return (Amp);
    else
      return (sigma);
  }

	
  void ABV(scalar r, scalar &A, scalar &B, scalar &V,
		  scalar &dAdr)
  {
    A = B = 1.0;
    dAdr = 0.0;

    V = Amp*exp(-r*r/(2.0*sigma*sigma));
  }

  scalar V(scalar r)
  {
    return (Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  scalar dVdr(scalar r)
  {
    return (-r/(sigma*sigma) * Amp*exp(-r*r/(2.0*sigma*sigma)));
  }

  scalar d2Adr2(scalar r)
  {
    return (0.0);
  }

  void Write (OutputSectionClass &outSection)
  {
    outSection.WriteVar ("Type", "Gaussian");
    outSection.WriteVar ("Amp", Amp);
    outSection.WriteVar ("sigma", sigma);
  }

  bool Read (InputSectionClass &inSection)
  {
    assert (inSection.ReadVar ("Amp", Amp));
    assert (inSection.ReadVar ("sigma", sigma));
  }

  PH_Gaussian(double Amp_, double sigma_)
  {
    Amp = Amp_;
    sigma = sigma_;
    // Do nothing for now
  }

  PH_Gaussian ()
  {
    Amp = 0.0;
    sigma = 1.0;
  }

  //void Write(char *FileName);
  //void Read (FILE *fin);
  
};



PseudoHamiltonian *ReadPH (InputBuffer &SectionBuf);
PseudoHamiltonian *ReadPH (InputSectionClass &inSection);
PseudoHamiltonian *Read_PH(char *FileName);

#endif
