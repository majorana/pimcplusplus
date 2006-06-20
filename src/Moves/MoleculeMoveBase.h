#ifndef MOLECULE_MOVE_BASE_H
#define MOLECULE_MOVE_BASE_H

//#include "MoveBase.h"
#include "StageClass.h"

//class MolMoveClass: public ParticleMoveClass{
class MolMoveClass: public LocalStageClass{
  public:
	int numAccepted, numMoves;
  // This will store an array of active particles for each molecule in the sim;
  // should eliminate the need to repeatedly asemble these arrays at each move
  Array<Array<int,1>,1> MolMembers;
  MolMoveClass(PathDataClass&, IO::IOSectionClass);
  dVec GetCOM(int slice, int mol);
  dVec TranslateMol(int slice, Array<int,1>& activePtcls, double epsilon);
  void RotateMol(int slice, Array<int,1>& activePtcls, dVec& axis, double theta);
  void RotateMol(int slice, Array<int,1>& activePtcls, double theta);
  void Read (IOSectionClass &in);
	void Accept();
};

dVec ArbitraryRotate(dVec axis,dVec coord, double phi);

#endif
