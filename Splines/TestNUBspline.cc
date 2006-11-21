#include "BsplineHelper.h"
#include "Grid.h"

void TestBasis()
{
  LinearGrid grid(0.0, 1.0, 11);
  NUBsplineBasis<LinearGrid> uniform;
  uniform.Init(&grid);

  FILE *fout = fopen ("uniformBasis.dat", "w");
  TinyVector<double,4> bfuncs;
  for (double x=0.0; x<=1.0; x+=0.0001) {
    uniform(x, bfuncs);
    fprintf (fout, "%22.16e %22.16e %22.16e %22.16e %22.16e\n",
	     x, bfuncs[0], bfuncs[1], bfuncs[2], bfuncs[3]);
  }
  fclose(fout);
}

void TestNUBasis()
{
  Array<double,1> points(11);
  points = 0.0, 1.0, 1.4, 1.8, 6.0, 7.0, 9.0, 10.0, 12.8, 13.0, 15.0;
  GeneralGrid grid;
  grid.Init(points);
  NUBsplineBasis<GeneralGrid> nonuni;
  nonuni.Init (&grid);

  FILE *fout = fopen ("NUBasis.dat", "w");
  TinyVector<double,4> bfuncs, dbfuncs, d2bfuncs;
  for (double x=0.0; x<=15.0; x+=0.001) {
//   for (int i=0; i<grid.NumPoints; i++) {
//     double x = grid(i);
    nonuni.Evaluate(x, bfuncs, dbfuncs, d2bfuncs);
    fprintf (fout, "%22.16e %22.16e %22.16e %22.16e %22.16e "
	     "%22.16e %22.16e %22.16e %22.16e "
	     "%22.16e %22.16e %22.16e %22.16e\n",
	     x, bfuncs[0], bfuncs[1], bfuncs[2], bfuncs[3],
	     dbfuncs[0], dbfuncs[1], dbfuncs[2], dbfuncs[3],
	     d2bfuncs[0], d2bfuncs[1], d2bfuncs[2], d2bfuncs[3]);
  }
  fclose(fout);
}

main()
{
  TestBasis();
  TestNUBasis();
}
