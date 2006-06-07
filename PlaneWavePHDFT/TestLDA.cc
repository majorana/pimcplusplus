#include "PlaneWavesMPI.h"
#include "../IO/IO.h"
#include "../MPI/Communication.h"

void TestInitCharge()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  in.OpenFile("Na_HF_NLPP.h5");
  Potential *pot = ReadPotential(in);
  in.CloseFile();

  int numBands  = 10;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  Array<Vec3,1> rions(16);
  rions(0)  = Vec3(  6.7701,   -0.5241,   -0.7601);
  rions(1)  = Vec3( -0.2356,  -11.1819,    0.1041);
  rions(2)  = Vec3(  0.7953,   -3.1228,   -3.3658);
  rions(3)  = Vec3( -9.6992,    2.5853,   -8.7713);
  rions(4)  = Vec3( -1.0770,  -12.9317,  -10.9011);
  rions(5)  = Vec3( -8.4076,   -4.4247,   -6.0560);
  rions(6)  = Vec3(-11.7379,    0.7186,    6.3867);
  rions(7)  = Vec3(  0.7048,    6.5211,   -7.6078);
  rions(8)  = Vec3( -0.9636,   -1.0799,    3.2954);
  rions(9)  = Vec3( 13.1805,   -9.5798,   -8.3901);
  rions(10) = Vec3( -1.8638,   -4.0041,  -11.5047);
  rions(11) = Vec3(  0.4408,    7.2936,    3.7794);
  rions(12) = Vec3(  8.0325,   -6.3283,    4.1347);
  rions(13) = Vec3(  5.4652,    2.5650,    5.4588);
  rions(14) = Vec3( -7.7732,   -6.0467,    1.7018);
  rions(15) = Vec3(  5.5497,  -11.6650,    9.0145);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 4.0, *pot, true, true);
  system.SetIons(rions);
  IOSectionClass out;
  out.NewFile ("InitDensity.h5");
  out.WriteVar("InitRho", system.GetDensity());
  out.CloseFile();
}


void TestSolveLDA()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  in.OpenFile("Na_HF_NLPP.h5");
  Potential *pot = ReadPotential(in);
  in.CloseFile();

  int numBands  = 14;
  int numElecs  = 16;
  Vec3 box (26.56, 26.56, 26.56);
  Array<Vec3,1> rions(16);
  rions(0)  = Vec3(  6.7701,   -0.5241,   -0.7601);
  rions(1)  = Vec3( -0.2356,  -11.1819,    0.1041);
  rions(2)  = Vec3(  0.7953,   -3.1228,   -3.3658);
  rions(3)  = Vec3( -9.6992,    2.5853,   -8.7713);
  rions(4)  = Vec3( -1.0770,  -12.9317,  -10.9011);
  rions(5)  = Vec3( -8.4076,   -4.4247,   -6.0560);
  rions(6)  = Vec3(-11.7379,    0.7186,    6.3867);
  rions(7)  = Vec3(  0.7048,    6.5211,   -7.6078);
  rions(8)  = Vec3( -0.9636,   -1.0799,    3.2954);
  rions(9)  = Vec3( 13.1805,   -9.5798,   -8.3901);
  rions(10) = Vec3( -1.8638,   -4.0041,  -11.5047);
  rions(11) = Vec3(  0.4408,    7.2936,    3.7794);
  rions(12) = Vec3(  8.0325,   -6.3283,    4.1347);
  rions(13) = Vec3(  5.4652,    2.5650,    5.4588);
  rions(14) = Vec3( -7.7732,   -6.0467,    1.7018);
  rions(15) = Vec3(  5.5497,  -11.6650,    9.0145);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 4.0, *pot, true, true);
  system.SetIons(rions);
  system.SolveLDA();
}

void TestSolidLDA()
{
  CommunicatorClass bandComm, kComm;
  bandComm.SetWorld();
  Array<int,1> root(1);
  root(0) = 0;
  bandComm.Subset (root, kComm);

  IOSectionClass in;
  //in.OpenFile("NaLocalPH.h5");
  in.OpenFile("Na_HF_NLPP.h5");
  Potential *pot = ReadPotential(in);
  in.CloseFile();

  int numBands  = 5;
  int numElecs  = 2;
  Vec3 box (8.11, 8.11, 8.11);
  Array<Vec3,1> rions(2);
  rions(0)  = Vec3( 0.00, 0.00, 0.00); 
  rions(1)  = Vec3( 4.05, 4.05, 4.05);
  Vec3 Gprim = 2.0*M_PI*Vec3(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  Vec3 k = 0.25 * Gprim;
  cerr << "k = " << k << endl;
  MPISystemClass system (numBands, numElecs, bandComm, kComm, true, false);
  
  system.Setup (box, k, 8.0, *pot, true, true);
  system.SetIons(rions);
  system.SolveLDA();
}


#include "FermiSmear.h"
void 
TestSmear()
{
  FILE *fout = fopen ("Smear.dat", "w");
  double mu = 0.0;
  MethfesselPaxton smearer;
  smearer.SetOrder(2);
  smearer.SetWidth(0.5);
  for (double E=-3.0; E<=3.0; E+=0.01) 
    fprintf (fout, "%1.12e %1.12e %1.12e\n", 
	     E, smearer.D(E, mu), smearer.S(E, mu));
  fclose(fout);
}	     


main(int argc, char **argv)
{
  COMM::Init(argc, argv);
  // TestInitCharge();
  // TestSmear();
  TestSolveLDA();
  // TestSolidLDA();
}
