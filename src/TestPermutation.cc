#include "PIMCClass.h"


void TestPerm(PIMCClass &pimc)
{
  PathDataClass &PathData = pimc.PathData;
  PathClass &Path = PathData.Path;

  pimc.PermuteTable.ConstructCycleTable(0, 0, 8);
  PermuteTableClass &PermuteTable = pimc.PermuteTable;

  cerr << "HTable = \n" << PermuteTable.HTable << endl;
  cerr << "CycleTable = \n";
  for (int i=0; i<PermuteTable.NumEntries; i++) {
    CycleClass &perm = PermuteTable.PermTable(i);
    cerr << "Ncycles = " << perm.Ncycles << endl;
    cerr << "Cycle = [";
    for (int j=0; j<perm.Ncycles; j++)
      cerr << perm.CycleRep[j] << " ";
    cerr << "]\n";
    cerr << "P = " << perm.P << " C = " << perm.C << endl;
    cerr << endl;
  }

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
  TestPerm (pimc);
}
