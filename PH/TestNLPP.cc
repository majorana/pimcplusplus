#include "NLPPClass.h"

void 
TestRead()
{
  IOSectionClass in;
  NLPPClass nlpp;


//   assert (in.OpenFile("/home/esler/NLPP/Na/CASINO/Na_CASINO_NLPP.h5"));  
//   nlpp.Read (in);
//   nlpp.SetupProjectors (4.0, 16.0);

  assert (in.OpenFile("/home/esler/NLPP/N/CASINO/N.h5"));
  nlpp.Read (in);
  nlpp.SetupProjectors(10.0, 40.0);

}

main()
{
  TestRead();
}
