#include "ConjGrad.h"

main()
{
  Vec3 box(25.0, 25.0, 25.0);
  Hamiltonian H(box, 2.0, 1.0);
  ConjGrad CG(H);
  clock_t start, end;
  start = clock();
  for (int i=0; i<30; i++)
    CG.Iterate();
  end = clock();

  fprintf (stderr, "Time = %1.3f\n", 
	   (double)(end-start)/(double)CLOCKS_PER_SEC);

}
