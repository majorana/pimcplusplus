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
  //  myPairAction.ReadDavidSquarerFile("He4.95.dm");
  MirroredArrayClass<int> myArray(5,5);
  myArray.Set(0,0,1);
}
