#include "Common/IO/InputOutput.h"


void TestHDF5Output()
{
  IOSectionClass HDF5out;

  Array<double,1> v(5);
  Array<double,2> m(3,3);
  v = 1.0, 2.0, 3.0, 4.0, 5.0;
  m = 
    1.0, 2.0, 3.0,
    4.0, 5.0, 5.0,
    7.0, 8.0, 9.0;
  Array<double,1> mplus;
  mplus = 10.0, 11.0, 12.0;
  Array<string,3> types(2,4,1);
  types = "Electron", "Proton", "Kenon", "Bryon", "Nonon", "Brainon",
    "Stepon", "Offon";

  HDF5out.NewFile ("HDF5test2.h5");
  HDF5out.WriteVar ("v", v);
  for (int i=0; i<1; i++)
    HDF5out.AppendVar ("v", 6.0+i);
  HDF5out.WriteVar ("m", m);
  HDF5out.AppendVar ("m", mplus);
  HDF5out.WriteVar ("types", types);
  HDF5out.NewSection ("Action");
  HDF5out.WriteVar ("tau", 0.1);
  HDF5out.WriteVar ("tau2", 0.1);
  HDF5out.CloseSection();
  HDF5out.NewSection("Action");
  HDF5out.NewSection ("PairAction");
  HDF5out.WriteVar ("Types", types);
  HDF5out.WriteVar ("mass", 1.0);
  HDF5out.CloseSection();
  HDF5out.CloseSection();
  HDF5out.NewSection("Action");
  HDF5out.CloseSection();
  HDF5out.NewSection("Moves");
  HDF5out.WriteVar ("frequency", 100.0);
  HDF5out.CloseSection();
  HDF5out.NewSection("System", "System.h5");
  HDF5out.WriteVar ("pi", 3.14159);
  HDF5out.CloseSection();
  HDF5out.PrintTree();
  HDF5out.CloseFile();
  

  cerr << "Finished writing.\n";

  IOSectionClass HDF5in;
  IOSectionClass *sec;
  HDF5in.OpenFile ("HDF5test2.h5");
  //HDF5in.IncludeSection ("Pi", "Pi2.h5");
  Array<string,3> inTypes;
  HDF5in.ReadVar ("types", inTypes);
  cerr << "types = " << inTypes << endl;

  assert(HDF5in.OpenSection ("Action"));
  HDF5in.PrintTree();
  cerr << endl << endl;
  double tau;
  assert(HDF5in.ReadVar ("tau", tau));
  cerr << "tau = " << tau << endl;
  HDF5in.CloseSection();

  HDF5in.PrintTree();
  HDF5in.CloseFile();
}


main()
{
  TestHDF5Output();
}
