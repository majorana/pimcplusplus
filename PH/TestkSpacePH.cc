#include "kSpacePH.h"

class kSpaceTest
{
protected:
  Vec3 Box;
  Vec3 kBox;
  Array<Vec3,1> kVecs;
  Array<double,2> H;
  double kCut;
  Potential &PH;
public:
  void SetBox(double x, double y, double z);
  void SetupkVecs(double kcut);
  void CalcHamiltonian ();
  void Diagonalize();
  kSpaceTest (Potential &ph) : PH(ph) 
  {
    // Do nothing for now
  }
};

void Test (IOSectionClass &in)
{
  string PHname;
  assert (in.ReadVar ("PHfile", PHname));
  IOSectionClass PHin;
  assert (PHin.OpenFile(PHname));
  Potential *PH = ReadPotential(PHin);
  PHin.CloseFile();

  Array<double,1> box;
  assert(in.ReadVar ("Box", box));
  assert(box.size() == 3);
  
}

main()
{
  cerr << "Hello world!\n";

}
