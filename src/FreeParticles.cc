#include "Common/Blitz.h"
#include "Common/Random/Random.h"
#include "Common/IO/InputOutput.h"

#define NDIM 3

typedef TinyVector<int,NDIM> State;
typedef TinyVector<double,3> dVec;

class ParticleClass
{
protected:
  int NumParticles;
  dVec Box;
  dVec kPrim;
  Array<State,1> AllowedStates;
  Array<complex<double>,1> Rhok;
  double kCut;
  double Esum, E2sum;
  int BlockSize, NumBlocks;
  double lambda, beta;
  CommunicatorClass Communicator;
  RandomClass Random;
public: 
  Array<State,1> OccupiedStates;
  virtual double AcceptProb (int ptclNum, State &newState) = 0;
  // Returns sampling probablilty ratio oldProb/NewProb;
  virtual double Sample(int &ptclToMove, State &newState);
  void MCStep();
  void Observe();
  void DumpResults();
  void Read (IOSectionClass &in);
  virtual void FillStates();
  void Run();
  ParticleClass() :
    Random(Communicator)
  {
    // Do nothing for now
  }
};


class BoltzmannonClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
};

class BosonClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
};

class FermionClass : public ParticleClass
{
public:
  double AcceptProb (int ptclNum, State &newState);
  void FillStates();
  void FillFermiSea();
};


////////////////////////////////////
// ParticleClass Member Functions //
////////////////////////////////////
void ParticleClass::Read(IOSectionClass &in)
{
  Array<double,1> box;
  assert(in.ReadVar ("Box", box));
  for (int i=0; i<NDIM; i++) {
    Box[i] = box(i);
    kPrim[i] = 2.0*M_PI/box(i);
  }
  assert(in.ReadVar ("NumParticles", NumParticles));
  assert(in.ReadVar ("kCut", kCut));
  assert(in.ReadVar ("beta", beta));
  assert(in.ReadVar ("BlockSize", BlockSize));
  assert(in.ReadVar ("NumBlocks", NumBlocks));
  Esum = 0.0;
  E2sum = 0.0;
  FillStates();
}

void ParticleClass::FillStates()
{
  OccupiedStates.resize(NumParticles);
  for (int ptcl=0; ptcl<NumParticles; ptcl++)
    OccupiedStates(ptcl) = 0;
}

void ParticleClass::MCStep()
{
  State newState;
  int ptclToMove;
 
  double sampleRatio = Sample (ptclToMove, newState);

  double toAccept = AcceptProb (ptclToMove, newState);
  if (toAccept*sampleRatio>Random.Local()){
    OccupiedStates(ptclToMove)=newState;
  }
  else {
  }
}

double ParticleClass::Sample(int &ptclToMove, State &newState)
{
  ptclToMove=Random.LocalInt(OccupiedStates.size());
  int toChange;
  if (Random.Local()>0.5)
    toChange=1;
  else
    toChange=-1;
  int dimToChange=Random.LocalInt(3);
  newState=OccupiedStates(toChange);
  newState[dimToChange]+=toChange;
  return 1.0;
}

void ParticleClass::Run()
{

}


///////////////////////////////////
// FermionClass Member Functions //
///////////////////////////////////

void FermionClass::FillStates()
{


}



double FermionClass::AcceptProb (int ptclToMove, State &newState)
{
  // First, reject if new state is already occupied -- Pauli exclusion.
  for (int ptcl=0; ptcl<NumParticles; ptcl++)
    if (ptcl != ptclToMove)
      if (newState == OccupiedStates(ptcl))
	return 0.0;

  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  return (exp (-beta * (newEnergy-oldEnergy)));
}


////////////////////////////////////////
// Boltzmannon Class Member Functions //
////////////////////////////////////////
double BoltzmannonClass::AcceptProb (int ptclToMove, State &newState)
{
  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  double acceptProb = exp (-beta * (newEnergy-oldEnergy));
  return acceptProb;
}



/////////////////////////////////
// BosonClass Member Functions //
/////////////////////////////////

double BosonClass::AcceptProb (int ptclToMove, State &newState)
{
  int Nold = 0;
  int Nnew = 0;
  for (int i=0; i<NumParticles; i++) {
    if (OccupiedStates(i) == newState)
      Nnew++;
    if (OccupiedStates(i) == OccupiedStates(ptclToMove))
      Nold++;
  }

  double newEnergy, oldEnergy;
  dVec newk, oldk;
  for (int i=0; i<NDIM; i++) {
    newk[i] = newState[i]*kPrim[i];
    oldk[i] = OccupiedStates(ptclToMove)[i] * kPrim[i];
  }
  
  newEnergy = lambda * dot (newk, newk);
  oldEnergy = lambda * dot (oldk, oldk);

  double acceptProb = exp (-beta * (newEnergy-oldEnergy));

  // Boson statistics:
  acceptProb *= (double)(Nnew+1)/(double)Nold;
  return (acceptProb);
}


//////////
// Main //
//////////
main(int argc, char **argv)
{

  if (argc < 2) {
    cout << "Usage:\n";
    cout << "FreeParticles myfile.in\n"; 
  }
  else {
    ParticleClass *particle;
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    string particleType;
    assert (in.ReadVar ("Type", particleType));
    if (particleType == "FERMION")
      particle = new FermionClass;
    else if (particleType == "BOSON")
      particle = new BosonClass;
    else if (particleType == "BOLTZMANNON")
      particle = new BoltzmannonClass;
    else {
      cerr << "Unrecognized particle type " << particleType << ". Exitting.\n";
      exit(1);
    }
      
    particle->Read (in);
    particle->Run();
    delete particle;
  }
}

