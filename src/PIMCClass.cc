#include "PIMCClass.h"
#include "Moves/BisectionMoveClass.h"
#include "Moves/OpenBisectionMoveClass.h"
#include "Moves/MetaMoves.h"
#include "Moves/BlockMove.h"
#include "Moves/BisectionBlock.h"
#include "Observables/ObservableClass.h"
#include <sstream>

void PIMCClass::Read(IOSectionClass &in)
{
  // Read the parallelization strategy
  PathData.Read (in);
  cerr << "Finished PathData Read.\n";
  // Read in the system information and set up the path
  assert(in.OpenSection("System"));
  PathData.Path.Read(in);
  in.CloseSection();
  cerr << "Finished Path read.\n";
  // Read in the action information
  assert(in.OpenSection("Action"));
  PathData.Action.Read(in);
  cerr << "Finished Action read.\n";
  PathData.Actions.Read(in);
  cerr << "Finished Actions read.\n";
  in.CloseSection();
  // Read in the Observables
  assert(in.OpenSection("Observables"));
  ReadObservables(in);
  cerr << "Finished Observables Read.\n";
  in.CloseSection();


  // Read in the Moves
  assert(in.OpenSection("Moves"));
  ReadMoves(in);
  cerr << "Finished Moves Read.\n";
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




void PIMCClass::ReadObservables(IOSectionClass &in)
{
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot= myProc==0;
  if (iAmRoot) {
    string outFileBase;
    assert(in.ReadVar("OutFileBase",outFileBase));
    ostringstream cloneNum;
    cloneNum << PathData.GetCloneNum();
    OutFileName = 
      outFileBase+ "." + cloneNum.str() + ".h5";
    OutFile.NewFile(OutFileName);
    OutFile.NewSection("RunInfo");
    RunInfo.Write(OutFile);
    OutFile.CloseSection();
    OutFile.NewSection("System");
    WriteSystemInfo();
    OutFile.CloseSection(); // "System" 
    OutFile.NewSection ("Observables");
  }
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
      if (iAmRoot)
	OutFile.NewSection(theObserveName);
      tempObs = new PairCorrelationClass(PathData,OutFile);
    }
    else if (theObserveType=="nofr"){
      assert(in.ReadVar("Name",theObserveName));
      if (iAmRoot)
	OutFile.NewSection(theObserveName);
      tempObs = new nofrClass(PathData,OutFile);
    }
    else if (theObserveType=="Energy"){
      if (iAmRoot)
	OutFile.NewSection("Energies");
      tempObs = new EnergyClass(PathData,OutFile);
    }
    else if (theObserveType=="PathDump"){
      if (iAmRoot)
	OutFile.NewSection("PathDump");
      tempObs=new PathDumpClass(PathData,OutFile);
    }
    else if (theObserveType=="WindingNumber"){
      if (iAmRoot)
	OutFile.NewSection("WindingNumber");
      tempObs=new WindingNumberClass(PathData,OutFile);
    }
    else if (theObserveType=="StructureFactor"){
      if (iAmRoot)
	OutFile.NewSection("StructureFactor");
      tempObs=new StructureFactorClass(PathData,OutFile);
    }
    else if (theObserveType=="Sign"){
      if (iAmRoot)
	OutFile.NewSection("Sigh");
      tempObs=new WeightClass(PathData,OutFile);
    }
    else {
      cerr<<"We do not recognize the observable "<<theObserveType<<endl;
	abort();
    }
    tempObs->Read(in);
    Observables(counter)=tempObs;
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();//Observable
  }
  if (iAmRoot)
    OutFile.CloseSection(); // "Observables"
}




void PIMCClass::ReadMoves(IOSectionClass &in)
{

  int numOfMoves=in.CountSections("Move");
  Moves.resize(numOfMoves);
  int steps;
  int myProc=PathData.Path.Communicator.MyProc();
  bool iAmRoot = myProc == 0;
  
  if (iAmRoot)
    OutFile.NewSection("Moves");
  for (int counter=0;counter<numOfMoves;counter++){
    string moveName;
    in.OpenSection("Move",counter);
    string MoveType;
    assert(in.ReadVar("type",MoveType));
    if (iAmRoot)
      OutFile.NewSection(MoveType);
    if (MoveType=="Bisection") {
      moveName = "BisectionMove";
      Moves(counter)=new BisectionMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="OpenBisection") {
      moveName = "OpenBisectionMove";
      Moves(counter)=new OpenBisectionMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="ShiftMove") {
      Moves(counter)=new ShiftMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="PrintMove") {
      Moves(counter)=new PrintMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="CycleBlock") {
      moveName = "CycleBlockMove";
      Moves(counter)=new CycleBlockMoveClass(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="PermMove") {
      moveName= "PermMove";
      Moves(counter)=new PermMove(PathData, OutFile);
      Moves(counter)->Read(in);
    }
    else if (MoveType=="BisectionBlock"){
      moveName="BisectionBlock";
      Moves(counter)=new BisectionBlockClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="RefSlice"){
      moveName="RefSlice";
      Moves(counter)=new RefSliceMoveClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else if (MoveType=="Displace"){
      moveName="Displace";
      Moves(counter)=new DisplaceMoveClass(PathData,OutFile);
      Moves(counter)->Read(in); 
    }
    else {
      cerr<<"This type of move is not recognized: "<< MoveType <<endl;
      abort();
    }
    if (iAmRoot)
      OutFile.CloseSection();
    in.CloseSection();
  }
  if (iAmRoot)
    OutFile.CloseSection (); // "Moves"
  
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
