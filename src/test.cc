#include "Common.h"
#include "PathClass.h"
#include "IdenticalParticlesClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "MirroredArrayClass.h"




void setupAction(ActionClass &myActionClass,ArrayOfIdenticalParticlesClass &myIDParticles)
{
  myActionClass.PairActionVector.resize(1);
  myActionClass.PairActionVector(0).ReadDavidSquarerFile("../inputs/testcoul.dm");
  myActionClass.PairMatrix.resize(2,2);
  myActionClass.PairMatrix=0;
  myActionClass.tau=myActionClass.PairActionVector(0).tau;
  myActionClass.myIdenticalParticleArray=&myIDParticles;

     
}


void setupIDParticleArray(PathDataClass &myPathData)
{
  ArrayOfIdenticalParticlesClass &myIDParticles=myPathData.IdenticalParticleArray;
  double tau=0.1; //This better be the same as in the squarer file! UGly!
  int NumTimeSlices=32;
  ElectronsClass myElectrons;
  ProtonsClass myProtons;
  myElectrons.NumParticles=1;
  myProtons.NumParticles=1;
  myElectrons.lambda=0.5;
  myProtons.lambda=0;
   
  myElectrons.Path.resize(1,NumTimeSlices);
  myProtons.Path.resize(1,NumTimeSlices);
  setMode(BOTHMODE);
  dVec zeroVector=0;
  for  (int counter=0;counter<myProtons.NumParticles;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      myProtons.Path.SetPos(counter,counter2,zeroVector);
    }
  }
  double sigma=sqrt(2*0.5*tau);
  dVec electronVector;
  for (int counter=0;counter<myElectrons.NumParticles;counter++){
    for (int counter2=0;counter2<NumTimeSlices;counter2++){
      electronVector=GaussianRandomVec(sigma);
      myElectrons.Path.SetPos(counter,counter2,electronVector);
    }
  }

  myPathData.IdenticalParticleArray.resize(2);
  myPathData.IdenticalParticleArray.Set(0,myElectrons);
  myPathData.IdenticalParticleArray.Set(1,myProtons);
  myPathData.NumTimeSlices=NumTimeSlices;
}
  

void setupMove(BisectionMoveClass &myBisectionMove,ShiftMove &myShiftMove, PathDataClass &thePathData)
{

  Array<int,1> ActiveSpecies(1);
  ActiveSpecies = 0;
  myBisectionMove.PathData=&thePathData;
  myBisectionMove.SetActiveSpecies(ActiveSpecies);
  myBisectionMove.SetNumParticlesToMove(1);
  myBisectionMove.StartTimeSlice=0;
  myBisectionMove.NumLevels=3;
  myShiftMove.PathData=&thePathData;

}






  
   





int main()

{

  PathDataClass myPathData;
  ActionClass myActionClass;
  setupIDParticleArray(myPathData);
  setupAction(myActionClass,myPathData.IdenticalParticleArray);
  BisectionMoveClass myBisectionMove;
  ShiftMove myShiftMove;
  setupMove(myBisectionMove,myShiftMove,myPathData);
  
  for (int counter=0;counter<10000;counter++){
    for (int counter2=0;counter2<2;counter2++){
      myBisectionMove.makeMove();
    }
    myShiftMove.makeMove();
  }


}
  
