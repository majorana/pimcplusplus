#ifndef SPECIES_CLASS_H
#define SPECIES_CLASS_H

#include "Common.h"
#include "Common/IO/InputOutput.h"

typedef enum {FERMION, BOSON, BOLTZMANNON, ANYON} ParticleType;


/// This is an base class that holds all the information about
/// identical particles.  It may be specialized to hold specialized
/// information about particular types of particles.
class SpeciesClass
{
public:
  string Name, Type, NodeType;
  /// FirstPtcl and LastPtcl are inclusive
  int LastPtcl;
  int FirstPtcl;
  int NumDim;
  int NumParticles;
  TinyVector <bool,NDIM> DimensionActive;
  virtual bool Read(IOSectionClass &inSection);
  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  double Charge;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual ParticleType GetParticleType() = 0;
};

SpeciesClass* ReadSpecies(IOSectionClass &inSection);



///This is the inherited class that holds the information about
///the electrons. Eventually this will probably be turned into
///a fermion class.  Currently, no actual information about the electrons
///are actually included.
class FermionClass : public SpeciesClass
{
public:
  bool Read(IOSectionClass &inSection);
  FermionClass()  { 
    NumDim=NDIM;  
    DimensionActive=true;
  };
  ~FermionClass() { };
  ParticleType GetParticleType(){ return FERMION; }
};


///This is the inherited class that holds the information about
///the protons.  Currently no useful information about protons is contained
///here
class BosonClass : public SpeciesClass
{
public:
  ParticleType GetParticleType(){ return BOSON; }
  ///Just returns 0 until we do something more intelligble.
  BosonClass() 
  {
    NumDim=NDIM;
    DimensionActive=true;
  }
  ~BosonClass()
  {
  }
  
};

class BoltzmannonClass : public SpeciesClass
{
public:
  ParticleType GetParticleType(){ return BOLTZMANNON; }

  BoltzmannonClass()
  {
    NumDim=NDIM;  
    DimensionActive=true;
  }
  
  ~BoltzmannonClass()
  {
  }
  
};

#endif

