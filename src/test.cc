#include "Common.h"
#include "PathClass.h"
#include "IdenticalParticlesClass.h"
#include "ActionClass.h"
#include "PathDataClass.h"
#include "BisectionMoveClass.h"
#include "MirroredArrayClass.h"




void setupAction(ActionClass &myActionClass,PathDataClass &myPathData)


{
  myActionClass.PairActionVector.resize(1);
  myActionClass.PairActionVector(0).ReadDavidSquarerFile("../inputs/testcoul.dm");
  myActionClass.PairMatrix.resize(2,2);
  myActionClass.PairMatrix=0;
  myActionClass.tau=myActionClass.PairActionVector(0).tau;
  myActionClass.myIdenticalParticleArray=&(myPathData.IdenticalParticleArray);

     
}


void setupIDParticleArray(PathDataClass &myPathData)
{
  ArrayOfIdenticalParticlesClass &myIDParticles=myPathData.IdenticalParticleArray;
  double tau=0.1; //This better be the same as in the squarer file! UGly!
  int NumTimeSlices=32;
  ElectronsClass *myElectronptr = new ElectronsClass;
  ProtonsClass *myProtonptr = new ProtonsClass;
  ElectronsClass &myElectrons = *myElectronptr;
  ProtonsClass &myProtons = *myProtonptr;
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
  ActiveSpecies(0) = 0;
  myBisectionMove.PathData=&thePathData;
  cerr << "ActiveSpecies = " << ActiveSpecies << endl;
  myBisectionMove.SetActiveSpecies(ActiveSpecies);
  myBisectionMove.SetNumParticlesToMove(1);
  myBisectionMove.StartTimeSlice=0;
  myBisectionMove.NumLevels=3;
  myShiftMove.PathData=&thePathData;

}






  
   


#include <mpi.h>


int main(int argc, char **argv)

{
  MPI_Init(&argc, &argv);
  PathDataClass myPathData;
  ActionClass myActionClass;
  setupIDParticleArray(myPathData);
  setupAction(myActionClass,myPathData);
  BisectionMoveClass myBisectionMove;
  ShiftMove myShiftMove;
  setupMove(myBisectionMove,myShiftMove,myPathData);
  cerr<<"The size of the IdenticalParticleArray is ";
  cerr << (myBisectionMove.PathData)->IdenticalParticleArray.size()<<endl;
  cerr<<"What the action class thinks the size is: ";
  cerr<<  myActionClass.myIdenticalParticleArray->size()<<endl;
  for (int counter=0;counter<10000;counter++){
    for (int counter2=0;counter2<2;counter2++){
      cerr << "Doing step " << counter << endl;
      myBisectionMove.makeMove();
    }
    myShiftMove.makeMove();
  }

  MPI_Finalize();
}
  
