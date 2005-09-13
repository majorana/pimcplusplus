#include "IO.h"

using namespace IO;

void TestHDF5()
{
  IOSectionClass out;
  
  out.NewFile("IOtest.h5");
  Array<double,2> myVar(2,3), myVar2;
  myVar(0,0)=1.0; myVar(0,1)=2.0; myVar(0,2)=3.0;
  myVar(1,0)=4.0; myVar(1,1)=5.0; myVar(1,2)=6.0;

  string TestString = "Hello World";
  out.WriteVar("TestString", TestString);

  out.WriteVar("Test", myVar);
  double a = 3.0;
  out.WriteVar ("a", a);

  out.CloseFile();

  Array<double,1> x(3);
  x(0) = 7.0; x(1) = 8.0; x(2)=9.0;

  IOSectionClass in;
  assert(in.OpenFile ("IOtest.h5"));
  IOVarBase *test = in.GetVarPtr("Test");
  nilArraySection n;
  test->Write(x,0,Range::all());
  in.PrintTree();
  in.ReadVar("Test", myVar2);
  cerr << "myVar2 = " << myVar2 << endl;
  Array<double,2> RangeVar;
  test->Read(RangeVar, Range(0,1), Range(0,2,2));
  cerr << "RangeVar = " << RangeVar << endl;
  in.CloseFile();

}


void TestASCII()
{
  cerr << "\n\nASCII Test:\n\n";
  Array<double,2> myVar;
  Array<double,1> y(2);
  Array<double,1> x(3);
  x(0) = 7.0; x(1) = 8.0; x(2)=9.0;

  IOSectionClass in;
  assert(in.OpenFile ("IOtest.txt"));
  in.PrintTree();
  IOVarBase *test = in.GetVarPtr("Test");
  nilArraySection n;
  test->Write(x,0,Range::all());
  test->Read(y,1,Range(1,2));
  cerr << "y = " << y << endl;
  in.PrintTree();
  in.ReadVar("Test", myVar);
  test->Read(myVar);
  cerr << "myVar = " << myVar << endl;
  //  in.CloseFile();

}


main()
{
  TestHDF5();
  TestASCII();
}
