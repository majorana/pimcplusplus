#include "Common.h"
#include "PathClass.h"
#include "IdenticleParticleClass.h"
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


void setupIDParticleArray(PathData &myPathData){

ArrayOfIdenticalParticlesClass &myIDParticles)
{

  double tau=0.1 //This better be the same as in the squarer file! UGly!
  int numTimeSlices=32;
  ElectronsClass myElectrons;
  ProtonsClass myProtons;
  myElectrons.NumParticles=1;
  myProtons.NumParticles=1;
  myElectrons.lambda=0.5;
  myProtons.lambda=0;
   
  myElectrons.Path.resize(1,numTimeSlices);
  myProtons.Path.resize(1,numTimeSlices);
  setMode(OBSERVABLEMODE);
  dVec zeroVector=0;
  for  (int counter=0;counter<myProtons.numParticles;counter++){
    for (int counter2=0;counter<numTimeSlices;counter2++){
      myProtons.Path.setPos(counter,counter2,zeroVector);
    }
  }
  double sigma=sqrt(2*0.5*tau);
  dVec electronVector;
  for (int counter=0;counter<myElectrons.numParticles;counter++){
    for (int counter2=0;counter2<numParticles;counter2++){
      electronVector=GaussianRandomVector;
      myElectrons.Path.setPos(counter,counter2,electronVector);
    }
  }

  myPathData.IdenticalParticlesArray.IdenticalParticleArray.resize(2);
  myPathData.IdenticalParticlesArray.IdenticalParticleArray(0)=&myElectrons;
  myPathData.IdenticalParticlesArray.IdenticalParticleArray(1)=&myProtons;
  myPathData.NumTimeSlices=numTimeSlices;
}
  

void setupMove(BisectionMoveClass &myBisectionMove,ShiftMove &myShiftMove, PathData &thePathData)
{

  Array<int,1> ActiveSpecies=0;
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
  
