#ifndef SPECIES_CLASS_H
#define SPECIES_CLASS_H


#include "PermutationClass.h"

typedef enum {FERMION, BOSON, BOLTZMANNON, ANYON} ParticleType;


/// This is an base class that holds all the information about
/// identical particles.  It may be specialized to hold specialized
/// information about particular types of particles.
class SpeciesClass
{
public:

  /// FirstPtcl and LastPtcl are inclusive
  int LastPtcl;
  int FirstPtcl;

  int NumParticles;

  /// \$ \lambda \equiv \frac{\hbar^2}{2m} \$.  This is zero for a
  /// classical particle.
  double lambda;
  
  /// Returns the nodal action for fermions.  Returns 0 for bosons.
  virtual double NodeAction (int Ptcl, int LinkNum) = 0;
  virtual ParticleType GetParticleType()=0;
};




///This is the inherited class that holds the information about
///the electrons. Eventually this will probably be turned into
///a fermion class.  Currently, no actual information about the electrons
///are actually included.
class FermionClass : public SpeciesClass
{
public:
  ///When we make this work, this calculates the NodeActions
  double NodeAction (int Ptcl, int LinkNum);
  FermionClass()  { };
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
  double NodeAction (int Ptcl, int LinkNum)
    {
      return (0.0);
    }  
  BosonClass()
    {

    }
  ~BosonClass()
    {
    }
  
};

class BoltzmannonClass : public SpeciesClass
{
public:
  ParticleType GetParticleType(){ return BOLTZMANNON; }
  ///Just returns 0 until we do something more intelligble.
  double NodeAction (int Ptcl, int LinkNum)
    {
      return (0.0);
    }  
  BoltzmannonClass()
    {

    }
  ~BoltzmannonClass()
    {
    }
  
};

#endif

