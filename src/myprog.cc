#include "PathClass.h"
#include "IdenticleParticleClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"

int main()

{
  PathClass myPath;
  ElectronsClass myElectrons;
  ActionClass myAction;
  PathDataClass myPathData;
  BisectionMoveClass myBisectionMove;
  PairActionClass myPairAction;
  myPairAction.ReadDavidSquarerFile("He.4.95.dm");
}
