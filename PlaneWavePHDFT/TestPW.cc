#include "ConjGrad.h"

main()
{
  Vec3 box(15.0, 15.0, 15.0);
  Hamiltonian H(box, 4.0, 1.0);
  ConjGrad CG(H);
  for (int i=0; i<3000; i++)
    CG.Iterate();

}
