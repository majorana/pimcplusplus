#include "PIMCClass.h"
#include "MirroredClass.h"


main(int argc, char **argv)
{
  string version = VERSION;
  cerr << "pimc++ v. " << version << endl;

  //  MirroredClassTest();
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  if (argc < 2) {
    cout << "Usage:\n";
    cout << "pimc++ myfile.in\n"; 
 }
  else {
    ifstream infile;
    infile.open(argv[1]);

    infile.seekg(0,ios::end);
    int length=infile.tellg();
    infile.seekg(0,ios::beg);
    char *buffer=new char[length];
    infile.read(buffer,length);
    infile.close(); 
    string fileCopy(buffer);   
    delete buffer;
    cerr<<"My buffer is "<<fileCopy<<endl;
    IOSectionClass in;
    assert (in.OpenFile(argv[1]));
    PIMCClass PIMC;
    PIMC.Read(in,fileCopy);
    PIMC.Run();
  }
}
