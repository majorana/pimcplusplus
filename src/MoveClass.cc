#include "MoveClass.h"


void ParticleMoveClass::SetActiveSpecies (Array<int,1> ActSpecies)
{
  ActiveSpecies.resize(ActSpecies.size());
  ActiveSpecies = ActSpecies;
  ///This calculates the total number of particles
  TotalParticles = 0;
  int CurrentNumPtcls = 1234567;
  for (int i=0; i<ActSpecies.size(); i++) {
    CurrentNumPtcls = 
      PathData.SpeciesArray(ActSpecies(i)).NumParticles(); 
    TotalParticles += CurrentNumPtcls;
  }


  ///This proceeds to calculate the pieces necessary for
  ///the particle id mapping..I don't think this is necessary 
  ///any longer,  but I'm not sure. 

  /*

  //  cerr << "TotalParticles = " << TotalParticles << endl;
  MyParticles.resize(TotalParticles);
  ///This sets the total number of particles 
  TotalParticles = 0;
  for (int i=0; i<ActSpecies.size(); i++) {
    int CurrentNumPtcls = 
      PathData.SpeciesArray(ActSpecies(i)).NumParticles(); 
    //    cerr << "CurrentNumPtcls = " << CurrentNumPtcls << endl;
    for (int j=0; j<CurrentNumPtcls; j++) {
      MyParticles(j+TotalParticles)[0] = ActSpecies(i);
      MyParticles(j+TotalParticles)[1] = j;
    }
    TotalParticles += CurrentNumPtcls;
  }

  */
}


inline int RandInt (int Max) //Hopefully this didn't break anything
{
  double myRandNum;
  double myNum;
  myRandNum=sprng();
  myNum=(double)Max * myRandNum;
    
  return ((int)floor(myNum));
}


/// So do we still want to choose particles by dumping everything
/// into some mapping array from teh active particles and dealing 
// with it that way?
void ParticleMoveClass::ChooseParticles()
{
  for (int i=0; i<NumParticlesToMove; i++) { 
    bool Redundant;
    do {
      MyParticleIndices(i) = RandInt(TotalParticles);
      Redundant = false;
      for (int j=0; j<i; j++){
	  if (MyParticleIndices(i) == MyParticleIndices(j)){
	      Redundant = true;
	      break;
	  }
      }      
    } while (Redundant); 
  }
  for (int i=0; i<NumParticlesToMove; i++) 
    ActiveParticles(i) = MyParticles(MyParticleIndices(i));
}
