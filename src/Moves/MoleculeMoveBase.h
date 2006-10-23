#ifndef MOLECULE_MOVE_BASE_H
#define MOLECULE_MOVE_BASE_H

//#include "MoveBase.h"
#include "StageClass.h"

// Identifies whether moves should be attempted globally
// or one particle at a time (SINGLE)
enum MoveMode{GLOBAL, SINGLE, SEQUENTIAL};

//class MolMoveClass: public ParticleMoveClass{
class MolMoveClass: public LocalStageClass{
	int numMol, molIndex;
  public:
	MoveMode mode;
	int numAccepted, numMoves;
	// located in PathClass now
  // This will store an array of active particles for each molecule in the sim;
  // should eliminate the need to repeatedly asemble these arrays at each move
  //Array<Array<int,1>,1> MolMembers;
	Array<int, 1> MoveList;
  MolMoveClass(PathDataClass&, IO::IOSectionClass);
  dVec GetCOM(int slice, int mol);
  dVec TranslateMol(int slice, Array<int,1>& activePtcls, double epsilon);
	void MolMoveClass::TranslatePtcl(int slice, int ptcl, double Sigma);
  void RotateMol(int slice, Array<int,1>& activePtcls, dVec& axis, double theta);
  void RotateMol(int slice, Array<int,1>& activePtcls, double theta);
  void Read (IOSectionClass &in);
	void Accept();
	void Advance();
};

dVec ArbitraryRotate(dVec axis,dVec coord, double phi);

#endif
