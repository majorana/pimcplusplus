#include "PAFit.h"


void PairActionFitClass::WriteHeader (Rho &rho, IOSectionClass &outSection)
{
  outSection.NewSection ("Particle1");
  Particle1.Write(outSection);
  outSection.CloseSection();

  outSection.NewSection ("Particle2");
  Particle2.Write(outSection);
  outSection.CloseSection();

  outSection.WriteVar ("UsePBC", UsePBC);
  if (UsePBC)
    outSection.WriteVar ("UsePBC");

  outSection.NewSection("PH");
  PH->Write (outSection);
}

void PA2DFitClass::Write(IOSectionClass &outSection)
{

}
