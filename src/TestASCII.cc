#include <string.h>
#include <list>
#include "Common/IO/InputOutput.h"

using namespace std;

int  main()
{

  IOSectionClass ASCIIout;

  Array<double,1> v(5);
  Array<double,3> m(3,3,1);
  v = 1.0, 2.0, 3.0, 4.0, 5.0;
  m = 
    1.0, 2.0, 3.0,
    4.0, 5.0, 6.0,
    7.0, 8.0, 9.0;
  Array<double,2> mplus(3,1);
  mplus = 10.0, 11.0, 12.0;
  Array<string,3> types(2,4,1);
  Array<bool,1> bools(3);
  bools = true, false, true;
  types = "Electron", "Proton", "Kenon", "Bryon", "Nonon", "Brainon",
    "Stepon", "Offon";
  Array<string,2> typeAdd(4,1);
  typeAdd="on", "off", "in", "out";
  ASCIIout.NewFile ("ASCIItest2.txt");
  ASCIIout.WriteVar ("v", v);
  ASCIIout.WriteVar ("bools", bools);
  for (int i=0; i<1; i++)
    ASCIIout.AppendVar ("v", 6.0+i);
  ASCIIout.WriteVar ("m", m);
  ASCIIout.AppendVar ("m", mplus);
  ASCIIout.ReadVar("m", m);
  cerr << "m = " << m << endl;
  ASCIIout.WriteVar ("types", types);
  ASCIIout.AppendVar("types",typeAdd);
  ASCIIout.NewSection ("Action");
  ASCIIout.WriteVar ("tau", 0.1);
  ASCIIout.WriteVar ("tau2", 0.1);
  ASCIIout.CloseSection();
  ASCIIout.NewSection("Action");
  ASCIIout.NewSection ("PairAction");
  ASCIIout.WriteVar ("Types", types);
  ASCIIout.WriteVar ("mass", 1.0);
  ASCIIout.FlushFile();
  ASCIIout.CloseSection();
  ASCIIout.CloseSection();
  ASCIIout.NewSection("Action");
  ASCIIout.CloseSection();
  ASCIIout.NewSection("Moves");
  ASCIIout.WriteVar ("frequency", 100.0);
  ASCIIout.CloseSection();
  ASCIIout.NewSection("System", "System.h5");
  ASCIIout.WriteVar ("pi", 3.14159);
  ASCIIout.CloseSection();
  ASCIIout.PrintTree();
  ASCIIout.CloseFile();
  

  cerr << "Finished writing.\n";
 



// //   IOTreeASCIIClass IO;
// //   IO.OpenFile("inputFile", "Root", NULL);
// //   IO.PrintTree();
//   IOSectionClass inSection;
//   inSection.OpenFile ("inputFile");
//   inSection.PrintTree();
//   double tau;
//   inSection.ReadVar("tau",tau);
//   cerr << "tau = " << tau << endl;
//   Array<double,3> Positions;
//   inSection.OpenSection("System");
//   inSection.OpenSection("Particles");
//   inSection.OpenSection("Species");
//   inSection.ReadVar ("Positions", Positions);
//   cerr << "Positions = " << Positions << endl;
//   inSection.CloseFile();
}
