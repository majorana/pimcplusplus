#include "PathDataClass.h"

/// FIX ME!  BAD BAD SLOW HACK!!!
void PathDataClass::acceptMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice)
{
  Array<int,1> Ptcl(1);
  for (int i=0; i<ActiveParticles.size(); i++)
    {
      int Species = ActiveParticles(i)[0];
      Ptcl(0) = ActiveParticles(i)[1];
      IdenticalParticleArray(Species).Path.AcceptCopy(Ptcl, StartTimeSlice,EndTimeSlice);
    }

}

void PathDataClass::rejectMove(Array <ParticleID,1> ActiveParticles,int StartTimeSlice,int EndTimeSlice)
{
  Array<int,1> Ptcl(1);
  for (int i=0; i<ActiveParticles.size(); i++)
    {
      int Species = ActiveParticles(i)[0];
      Ptcl(0) = ActiveParticles(i)[1];
      IdenticalParticleArray(Species).Path.RejectCopy(Ptcl, StartTimeSlice,EndTimeSlice);
    }
}
