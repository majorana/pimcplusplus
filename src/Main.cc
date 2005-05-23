#include "PIMCClass.h"
#include "MirroredClass.h"


main(int argc, char **argv)
{
  string version = VERSION;
  cerr << "pimc++ v. " << version << endl;

  //  MirroredClassTest();
  COMM::Init(argc, argv);
  if (argc < 2) {
    cout << "Usage:\n";
    cout << "pimc++ myfile.in\n"; 
 }
  else {
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    PIMCClass PIMC;
    PIMC.Read(in);
    PIMC.Run();
  }
}
