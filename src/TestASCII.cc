#include <string.h>
#include <list>
#include "Common/IO/InputOutput.h"

using namespace std;

int  main()
{
//   InputTreeASCIIClass IO;
//   IO.OpenFile("inputFile", "Root", NULL);
//   IO.PrintTree();
  InputSectionClass inSection;
  inSection.OpenFile ("inputFile");
  inSection.PrintTree();
  double tau;
  inSection.ReadVar("tau",tau);
  cerr << "tau = " << tau << endl;
  Array<double,3> Positions;
  inSection.OpenSection("System");
  inSection.OpenSection("Particles");
  inSection.OpenSection("Species");
  inSection.ReadVar ("Positions", Positions);
  cerr << "Positions = " << Positions << endl;

}
