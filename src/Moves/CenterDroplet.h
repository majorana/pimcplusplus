#ifndef CENTERDROPLET_MOVE_H
#define CENTERDROPLET_MOVE_H

#include "MoveBase.h"


class CenterDropletClass : public ParticleMoveClass
{
 public:
  int numAccepted,numMoves;
  void MakeMove();
  int Species;

  inline void Read(IOSectionClass &moveInput)
    {
      string typeCheck;
      assert(moveInput.ReadVar("type",typeCheck));
      assert(typeCheck=="CenterDroplet");
      assert(moveInput.ReadVar("name",Name));
      string speciesString;
      assert(moveInput.ReadVar("species",speciesString));
      Species=PathData.Path.SpeciesNum(speciesString);

    }
  CenterDropletClass(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection),numAccepted(0),numMoves(0)
    {
      ActiveParticles.resize(PathData.Path.NumParticles());
      for (int ptcl=0;ptcl<PathData.Path.NumParticles();ptcl++)
	ActiveParticles(ptcl)=ptcl;
      /* do nothing for now */
    }
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};



#endif
