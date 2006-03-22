#include "Grid.h"
#include "../IO/IO.h"

void TestOptimal()
{
  OptimalGrid rgrid(600, 15.0, 0.0039);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Optimal ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
}


void TestOptimal2()
{
  OptimalGrid2 rgrid(1.0e-6, 15, 10.0, 300);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Optimal2 ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
}


void TestCluster()
{
  ClusterGrid rgrid(1.0e-6, 15, 1.0/10.0, 300);
  double delta = rgrid(1) - rgrid(0);

  for (int i=0; i<rgrid.NumPoints; i++) {
    double r = rgrid(i) + 0.1*delta;
    int j = rgrid.ReverseMap(r);
    if (i != j)
      cerr << "Error in Cluster ReverseMap!\n"
	   << " i = " << i << "j = " << j << endl;
  }
  for (int i=0; i<rgrid.NumPoints; i++)
    fprintf (stdout, "%1.16e\n", rgrid(i));
}


void TestGeneral()
{
  IOSectionClass in;
  in.OpenFile("Na_HF_NLP.h5");
  in.OpenSection("Grid");
  Grid *grid = ReadGrid(in);
  int lo = grid->ReverseMap(3.2);
  cerr << "Grid(lo) = " << (*grid)(lo) << endl;
  cerr << "Grid(hi) = " << (*grid)(lo+1) << endl;
  in.CloseFile();
}


main()
{
  TestGeneral();
//   TestOptimal();
//   TestOptimal2();
//   TestCluster();
}
