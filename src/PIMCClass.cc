#include "PIMCClass.h"
#include "BisectionMoveClass.h"
#include "OpenBisectionMoveClass.h"
#include "MetaMoves.h"
#include "BlockMove.h"
#include "Observables/ObservableClass.h"


void PIMCClass::Read(IOSectionClass &in)
{
  // Read in the system information and set up the path
  assert(in.OpenSection("System"));
  PathData.Path.Read(in);
  in.CloseSection();

  // Read in the action information
  assert(in.OpenSection("Action"));
  PathData.Action.Read(in);
  in.CloseSection();

  // Read in the Moves
  assert(in.OpenSection("Moves"));
  ReadMoves(in);
  in.CloseSection();

  // Read in the Observables
  assert(in.OpenSection("Observables"));
  ReadObservables(in);
  in.CloseSection();
  
  // Read in the Algorithm
  assert(in.OpenSection("Algorithm"));
  ReadAlgorithm(in);
  in.CloseSection();

  // Read in the Permuation info
//   assert(in.OpenSection("Permutations"));
//   ForwPermuteTable.Read(in);
//   RevPermuteTable.Read(in);
//   in.CloseSection();
}


void PIMCClass::ReadMoves(IOSectionClass &in)
{

  int numOfMoves=in.CountSections("Move");
  Moves.resize(numOfMoves);
  int steps;
  for (int counter=0;counter<numOfMoves;counter++){
    in.OpenSection("Move",counter);
    string MoveType;
    assert(in.ReadVar("type",MoveType));
    if (MoveType=="Bisection")
      Moves(counter)=new BisectionMoveClass(PathData);
    else if (MoveType=="OpenBisection")
      Moves(counter)=new OpenBisectionMoveClass(PathData);
    else if (MoveType=="ShiftMove")
      Moves(counter)=new ShiftMoveClass(PathData);
    else if (MoveType=="PrintMove")
      Moves(counter)=new PrintMoveClass(PathData);
    else if (MoveType=="CycleBlock")
      Moves(counter)=new CycleBlockMoveClass(PathData);
    else if (MoveType=="PermMove")
      Moves(counter)=new PermMove(PathData);
    else {
      cerr<<"This type of move is not recognized: "<< MoveType <<endl;
      abort();
    }
    Moves(counter)->Read(in);
    in.CloseSection();
  }

}

void PIMCClass::ReadObservables(IOSectionClass &in)
{
  assert(in.ReadVar("OutFile",OutFileName));
  OutFile.NewFile(OutFileName);
  OutFile.NewSection("RunInfo");
  RunInfo.Write(OutFile);
  OutFile.CloseSection();
  OutFile.NewSection("System");
  WriteSystemInfo();
  OutFile.CloseSection(); // "System" 
  int numOfObservables=in.CountSections("Observable");
  Observables.resize(numOfObservables);
  for (int counter=0;counter<numOfObservables;counter++){
    in.OpenSection("Observable",counter);
    string theObserveType;
    string theObserveName;
    assert(in.ReadVar("type",theObserveType));
    ObservableClass* tempObs;
    if (theObserveType=="PairCorrelation"){
      assert(in.ReadVar("Name",theObserveName));
      OutFile.NewSection(theObserveName);
      tempObs = new PairCorrelationClass(PathData,OutFile);
    }
    else if (theObserveType=="nofr"){
      assert(in.ReadVar("Name",theObserveName));
      OutFile.NewSection(theObserveName);
      tempObs = new nofrClass(PathData,OutFile);
    }
    else if (theObserveType=="Energy"){
      OutFile.NewSection("Energies");
      tempObs = new TotalEnergyClass(PathData,OutFile);
    }
    else if (theObserveType=="PathDump"){
      OutFile.NewSection("PathDump");
      tempObs=new PathDumpClass(PathData,OutFile);
    }
    else if (theObserveType=="WindingNumber"){
      OutFile.NewSection("WindingNumber");
      tempObs=new WindingNumberClass(PathData,OutFile);
    }
    else {
      cerr<<"We do not recognize the observable "<<theObserveType<<endl;
      abort();
    }
    tempObs->Read(in);
    Observables(counter)=tempObs;
    OutFile.CloseSection();
    in.CloseSection();//Observable
  }
  

}

void PIMCClass::ReadAlgorithm(IOSectionClass &in)
{
  Algorithm.Read(in,1);
}


void PIMCClass::Run()
{
  Algorithm.DoEvent();
  //  Array<MoveClass*,1> Moves;
  cerr<<"hello"<<endl;
  for (int counter=0;counter<Moves.size();counter++){
    cout<<"My name is "<<((MoveClass*)Moves(counter))->Name<<endl;
    cout<<"My acceptance ratio is "<<((MoveClass*)Moves(counter))->AcceptanceRatio()<<endl;
  }
  
}

void PIMCClass::WriteSystemInfo()
{
  dVec box = PathData.Path.GetBox();
  Array<double,1> boxArray(3);
  boxArray(0) = box[0];   boxArray(1) = box[1];   boxArray(2) = box[2];
  OutFile.WriteVar ("Box", boxArray);
  OutFile.WriteVar("tau",PathData.Action.tau);
  OutFile.WriteVar("NumTimeSlices",PathData.Path.TotalNumSlices);
  for (int speciesIndex=0; speciesIndex < PathData.Path.NumSpecies(); 
       speciesIndex++) {
    SpeciesClass &species = PathData.Path.Species(speciesIndex);
    OutFile.NewSection("Species");
    OutFile.WriteVar ("Name", species.Name);
    OutFile.WriteVar ("NumParticles", species.NumParticles);
    OutFile.WriteVar ("lambda", species.lambda);
    ParticleType type = species.GetParticleType();
    if (type == FERMION)
      OutFile.WriteVar ("ParticleType", "Fermion");
    if (type == BOSON)
      OutFile.WriteVar ("ParticleType", "Boson");
    if (type == BOLTZMANNON)
      OutFile.WriteVar ("ParticleType", "Boltzmannon");
    OutFile.CloseSection(); //"Species"
  }
}
