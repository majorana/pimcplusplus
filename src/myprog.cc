#include "Common.h"
#include "PathClass.h"
#include "IdenticleParticleClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "MirroredArrayClass.h"

int main()

{
  
  PathClass myPath;
  ElectronsClass myElectrons;
  ActionClass myAction;
  PathDataClass myPathData;
  BisectionMoveClass myBisectionMove;
  PairActionClass myPairAction;
  CommClass myCommunicator;
  //  myPairAction.ReadDavidSquarerFile("He4.95.dm");
  MirroredArrayClass<int> myArray(5,5);
  for (int counter=0;counter<5;counter++){
    for (int counter2=0;counter2<5;counter2++){
      myArray.Set(counter,counter2,counter*10+counter2);
    }
  }
  myArray.Print();
  myArray.shiftData(-3,myCommunicator);
  cout<<endl<<endl;
  myArray.Print();


}
