#include "MoveClass.h"

void SetActiveSpecies (Array<int,1> ActSpecies)
{
  ActiveSpecies.resize(ActSpecies.size());
  ActiveSpecies = ActSpecies;

  TotalParticles = 0;
  for (int i=0; i<ActSpecies.size(); i++) {
    int CurrentNumPtcls = 
      PathData->IdenticalParticleArray(ActSpecies(i)).NumParticles; 

    TotalParticles += CurrentNumPtcls;
  }

  MyParticles.resize(TotalParticles);
  
  TotalParticles = 0;
  for (int i=0; i<ActSpecies.size(); i++) {
    int CurrentNumPtcls = 
      PathData->IdenticalParticleArray(ActSpecies(i)).NumParticles; 

    for (int j=0; j<CurrentNumPtcls; j++) {
      MyParticles(j+TotalParticles)[0] = ActSpecies(i);
      MyParticles(j+TotalParticles)[1] = j;
    }
    TotalParticles += CurrentNumPtcls;
  }
}


inline int RandInt (int Max)
{
  return (floor((double)Max*sprng()));
}



void MoveClass::ChooseParticles(Array<ParticleID,1> &Particles)
{
  for (int i=0; i<NumParticlesToMove; i++) { 
    bool Redundant;
    do {
      MyParticleIndices(i) = RandInt(TotalParticles);
      Redundant = false;
      for (int j=0; j<i; j++)
	{
	  if (MyParticleIndices(i) == MyParticleIndices(j))
	    {
	      Redundant = true;
	      break;
	    }
	}      
    } while (Redundant); 
  }
  for (int i=0; i<NumParticlesToMove; i++) 
    ActiveParticles(i) = MyParticles(MyParticleIndices(i));
}
