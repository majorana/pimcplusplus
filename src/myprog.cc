#include "Common.h"
#include "PathClass.h"
#include "SpeciesClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "Moves/BisectionMoveClass.h"
#include "MirroredClass.h"

using namespace std;
int main()

{
  
  PathClass myPath;
  ElectronsClass myElectrons;
  ActionClass myAction;
  PathDataClass myPathData;
  BisectionMoveClass myBisectionMove;
  PairActionClass myPairAction;
  CommunicatorClass myCommunicator;
  myPairAction.ReadDavidSquarerFile("../inputs/testcoul.dm");

  LinearGrid grid(0.0, 3.0, 31);
  for (int i=0; i<grid.NumPoints; i++)
    {
      double r = grid(i);
      double s,q,z;
      // Calc diagonal action
      s = 0.0;
      q = r;
      z = 0.0;
      double UD = myPairAction.calcUsqz(s,q,z,6);
      // Calc far-off-diagonal action
      s = 2.0*r;
      q = r;
      z = 0.0;
      double UOD = myPairAction.calcUsqz(s,q,z,6);
      fprintf (stderr, "%4.1f %8.5f %8.5f\n", r, UD, UOD);
    }

  /*  MirroredArrayClass<int> myArray(5,5);
  for (int counter=0;counter<5;counter++){
    for (int counter2=0;counter2<5;counter2++){
      myArray.Set(counter,counter2,counter*10+counter2);
    }
  }
  myArray.Print();
  myArray.shiftData(-3,myCommunicator);
  cout<<endl<<endl;
  myArray.Print();
  */

}
