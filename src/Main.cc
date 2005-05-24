#include "PIMCClass.h"
#include "MirroredClass.h"

///  \mainpage An overview of PIMC++ 
///
/// \section intro_sec Introduction
///
/// Understanding PIMC++ can be broken down into understanding four
/// main pieces -- 
/// the moves, the observables, and the actions and how they rest on
/// top of the representation of the path.  
///
///  \section path_sec  The Path
/// There is a file called PathClass.cc that contains all the
/// information relevant to the path. PathClass stores information
/// concerning --the location of all the particles at all places in
/// imaginary time
/// --the species information
/// --the box information
///
/// \subsection species_subsection The Species
/// Every particle is assigned a specific species.  Each species has
/// certain information associated with it. This includes
///  
/// \section moves_sec The Moves
/// There are two ways that a move can be built in PIMC++
///

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
