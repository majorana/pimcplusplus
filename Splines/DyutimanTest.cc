#include "MyTricubicSpline.h"

void DyutimanTest()
{
  // Read data.
  FILE *fin;
  assert((fin = fopen ("orbitals.dat", "r")) != NULL);

  Array<double, 3> psi(71,45,19);
  
  for (int k=0; k<19; k++)
    for (int j=0; j<45; j++)
      for (int i=0; i<71; i++)
	fscanf (fin, " %lf ", &(psi(i,j,k)));
  fclose(fin);

  LinearGrid xgrid (0.0, 1.0, 71);
  LinearGrid ygrid (0.0, 1.0, 45);
  LinearGrid zgrid (0.0, 1.0, 19);

  MyTricubicSpline spline;
  spline.Init (&xgrid, &ygrid, &zgrid, psi);

  // Now do a test in the x direction.
  LinearGrid xtest(0.0, 1.0, 1000);
  double x, y, z;
  y = 0.444444; z = 0.31579;
  for (int i=0; i<1000; i++) {
    x = xtest(i);
    double psi = spline (x, y, z);
    double laplacian = spline.Laplacian(x,y,z);
    fprintf (stderr, "%1.12e %1.12e %1.12e\n", x, psi, laplacian);
  }
}


main()
{
  DyutimanTest();

}
