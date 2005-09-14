#include "IO.h"

using namespace IO;


void TestDoubleWrite(IOSectionClass &io)
{
  io.NewSection("Double");

  double x0 = 2.0;
  io.WriteVar("x0", x0);

  Array<double,1> x1(5);
  x1 =  1.0, 2.0, 3.0, 4.0, 5.0 ;
  io.WriteVar("x1", x1);

  Array<double,2> x2(5,3);
  x2 =  
    1.0,   2.0,  3.0, 
    4.0,   5.0,  5.0,
    7.0,   8.0,  9.0,
    10.0, 11.0, 12.0,
    23.0, 14.0, 15.0;
  io.WriteVar("x2", x2);
  
  IOVarBase* ioVar = io.GetVarPtr("x2");
  ioVar->Write(x1,Range::all(), 0);

  io.CloseSection();

  io.FlushFile();
}


void TestDoubleRead(IOSectionClass &io)
{
  io.OpenSection("Double");

  double x0;
  assert(io.ReadVar("x0", x0));
  assert (fabs(x0-2.0) < 1.0e-15);

  Array<double,1> x1;
  assert(io.ReadVar("x1", x1));
  for (int i=0; i<5; i++)
    assert (fabs(x1(i) - (i+1.0)) < 1.0e-15);

  Array<double,2> x2, x2test(5,3);
  x2test =  
    1.0,  2.0,  3.0, 
    2.0,  5.0,  5.0,
    3.0,  8.0,  9.0,
    4.0, 11.0, 12.0,
    5.0, 14.0, 15.0;
  
  assert(io.ReadVar("x2", x2));
  cerr << " x2 = " << x2 << endl;
  cerr << " x2test = " << x2test << endl;
  for (int i=0; i<5; i++)
    for (int j=0; j<3; j++) {
      assert (fabs(x2(i,j)-x2test(i,j)) < 1.0e-15);
    }

  Array<double,1> x2Range;
  IOVarBase* ioVar = io.GetVarPtr("x2");
  ioVar->Read(x2Range,Range::all(), 1);
  
  io.CloseSection();
}




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
//   TestHDF5();
//   TestASCII();
  IOSectionClass ioHDF5;
  ioHDF5.NewFile ("TestHDF5.h5");
  TestDoubleWrite(ioHDF5);
  TestDoubleRead(ioHDF5);

  IOSectionClass ioASCII;
  ioASCII.NewFile ("TestASCII.txt");
  TestDoubleWrite(ioASCII);
  TestDoubleRead(ioASCII);

}
