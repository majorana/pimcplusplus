#include "PathDataClass.h"

/// FIX ME!  BAD BAD SLOW HACK!!!
void PathDataClass::AcceptMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice)
{
  //  Array<int> NumPerSpecies;
  //  int speciesSize=PathDataClass.IdenticalParticleArrray.size();
  //  NumPerSpecies.resize(speciesSize); 
  //  NumPerSpecies=0;
  //  Array<Array<int,1>,1> ParticleSpeciesArray;
  //  ParticleSpeciesArray.resize(speciesSize);
  //  for (counter=0;counter<
  //  for (int i=0;i<

  for (int i=0; i<ActiveParticles.size(); i++)
    {
      int Species = ActiveParticles(i)[0];
      int Ptcl = ActiveParticles(i)[1];
      SpeciesArray(Species).Path.AcceptCopy(Ptcl, StartTimeSlice,
						      EndTimeSlice);
    }

}

void PathDataClass::RejectMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice)
{
  int Ptcl;
  for (int i=0; i<ActiveParticles.size(); i++)
    {
      int Species = ActiveParticles(i)[0];
      Ptcl = ActiveParticles(i)[1];
      SpeciesArray(Species).Path.RejectCopy(Ptcl, StartTimeSlice,
						      EndTimeSlice);
    }
}
