#include "PIMCClass.h"
#include "MirroredClass.h"


main(int argc, char **argv)
{
  //  MirroredClassTest();
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
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
