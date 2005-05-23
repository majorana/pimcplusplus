#include "PIMCClass.h"
#include "MirroredClass.h"


main(int argc, char **argv)
{
  //  MirroredClassTest();
  COMM::Init(argc, argv);
  string version = VERSION;
  perr << "pimc++ v. " << version << endl;

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
