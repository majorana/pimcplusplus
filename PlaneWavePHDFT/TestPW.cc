#include "ConjGrad.h"

main()
{
  Vec3 box(20.0, 20.0, 20.0);
  Hamiltonian H(box, 2.0, 1.0);
  ConjGrad CG(H);
  for (int i=0; i<30; i++)
    CG.Iterate();

}
