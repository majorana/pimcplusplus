#include "PIMCClass.h"


 void TestPerm(PIMCClass &pimc)
{
//   PathDataClass &PathData = pimc.PathData;
//   PathClass &Path = PathData.Path;

//   pimc.ForwPermuteTable.ConstructCycleTable(0, 0, 8);
//   PermuteTableClass &ForwPermuteTable = pimc.ForwPermuteTable;

//   cerr << "HTable = \n" << ForwPermuteTable.HTable << endl;
//   cerr << "CycleTable = \n";
//   for (int i=0; i<ForwPermuteTable.NumEntries; i++) {
//     CycleClass &cycle = ForwPermuteTable.CycleTable(i);
//     cerr << "Length = " << cycle.Length << endl;
//     cerr << "Cycle = [";
//     for (int j=0; j<cycle.Length; j++)
//       cerr << cycle.CycleRep[j] << " ";
//     cerr << "]\n";
//     cerr << "P = " << cycle.P << " C = " << cycle.C << endl;
//     cerr << endl;
//   }
//   //  cerr<<"MAde it here"<<endl;
  
//   // double forwardProb=pimc.ForwPermuteTable.AttemptPermutation();
//   pimc.ForwPermuteTable.CurrentCycle=pimc.ForwPermuteTable.CycleTable(1);
//   Array<int,1> ActPtcls(2);
//   ActPtcls(0)=0;
//   ActPtcls(1)=1;

//   pimc.Observables(0)->WriteBlock();
//   pimc.PathData.MoveJoin(3);
//   pimc.ForwPermuteTable.CurrentCycle.Apply(PathData.Path,0,3);
//   pimc.PathData.AcceptMove(2,4,ActPtcls);
//   pimc.Observables(0)->WriteBlock();
//   pimc.OutFile.CloseFile();
  
//   //  cerr<<"calculated forward prob"<<endl;
//   double reverseProb=
//     pimc.RevPermuteTable.CalcReverseProb(ForwPermuteTable);
//   //  cerr<<forwardProb<<" "<<reverseProb << endl;

//   int N = ForwPermuteTable.NumEntries;
//   Array<int,1> hits(N);
//   hits  = 0;
//   for (int i=0; i<10000000; i++) {
//     double xi = pimc.PathData.Path.Random.Local();
//     hits(ForwPermuteTable.FindEntry(xi))++;
//   }
  
//   for (int i=0; i<N; i++) {
//     double expProb = 
//       ForwPermuteTable.CycleTable(i).P * ForwPermuteTable.NormInv;
//     double obsProb = (double)hits(i) * 1.0e-6;
//     cerr << "Expected prob: "<<expProb<<endl;
//     cerr << "Observed prob: "<<obsProb<<endl;
//   }
    

}


main(int argc, char **argv)
{
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  PIMCClass pimc;
  IOSectionClass in;
  in.OpenFile ("TestPerm.in");
  pimc.Read(in);
  pimc.Run();
  TestPerm (pimc);
  
}
