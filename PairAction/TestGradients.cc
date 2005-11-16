#include "PAFit.h"
#include "../IO/FileExpand.h"

void TestGradients (IOSectionClass &in)
{
  string PAfilename;
  IOSectionClass PAIO;
  assert (in.ReadVar("PAfilename", PAfilename));
  PAfilename = ExpandFileName(PAfilename);
  assert(PAIO.OpenFile(PAfilename));
  double tau;
  assert (in.ReadVar("tau", tau));
  PairActionFitClass &PA = *(ReadPAFit(PAIO, tau, 1));
 
  PAIO.CloseFile();
}


main(int argc, char **argv)
{
  if (argc < 2) {
    cerr << "Usage:\n"
	 << "  TestGradients myfile.in\n";
    exit(1);
  }
  else {
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    TestGradients(in);
  }
}
