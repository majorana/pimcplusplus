#include "InputOutputHDF5.h"

void TestHDF5Output()
{
  OutputSectionHDF5Class HDF5out;

  Array<double,1> v(5);
  Array<double,2> m(3,3);
  v = 1.0, 2.0, 3.0, 4.0, 5.0;
  m = 
    1.0, 2.0, 3.0,
    4.0, 5.0, 5.0,
    7.0, 8.0, 9.0;

  HDF5out.OpenFile ("HDF5test.h5");
  HDF5out.WriteVar ("v", v);
  HDF5out.WriteVar ("m", m);
  HDF5out.OpenSection ("Action");
  HDF5out.WriteVar ("tau", 0.1);
  HDF5out.CloseSection();
  HDF5out.OpenSection("Action");
  HDF5out.OpenSection ("PairAction");
  HDF5out.WriteVar ("mass", 1.0);
  HDF5out.CloseSection();
  HDF5out.CloseSection();
  HDF5out.OpenSection("Action");
  HDF5out.CloseSection();
  HDF5out.OpenSection("Moves");
  HDF5out.WriteVar ("frequency", 100.0);
  HDF5out.CloseSection();
  HDF5out.CloseFile();

  InputSectionHDF5Class HDF5in;
  HDF5in.OpenFile ("HDF5test.h5",NULL);
  
  HDF5in.CloseFile();
}


main()
{
  TestHDF5Output();
}
