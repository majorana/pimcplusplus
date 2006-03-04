#ifndef CENTER_OF_MASS_MOVE_H
#define CENTER_OF_MASS_MOVE_H

#include "MoveBase.h"

class CenterOfMassMoveClass : public ParticleMoveClass
{
 public:
  void MakeMove();
  void Read(IOSectionClass &moveInput);
  CenterOfMassMoveClass(PathDataClass &myPathData,IOSectionClass outSection) : 
    ParticleMoveClass (myPathData,outSection)
    {

    }
  dVec original_center_of_mass;
  bool firstTime;
  //  double AcceptanceRatio(int numAccepted,int numMoves);
  void WriteRatio()
    {
      //do nothing for now
    };



};




#endif
