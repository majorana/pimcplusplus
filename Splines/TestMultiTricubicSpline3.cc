#include "MyTricubicSpline.h"
#include "MultiTricubicSpline3.h"

void ValTest()
{
  int N = 10;
  const int numSplines = 4;
  Array<MyTricubicSpline,1> MySplines(numSplines);
  MultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<double,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++) 
	  initData(ix,iy,iz,n) = drand48();

  cerr << "initData filled.\n";
  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp;
    temp.reference (initData(Range::all(),Range::all(),Range::all(),n));
    MySplines(n).Init (&xGrid, &yGrid, &zGrid, temp);
  }
  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  FILE *f1 = fopen ("tri1.dat", "w");
  FILE *f2 = fopen ("tri2.dat", "w");

  Array<double,1> vals(numSplines);
  for (double u=0.0; u<=1.0; u+=0.001) {
    double x = xGrid.Start + u*(xGrid.End-xGrid.Start);
    double y = xGrid.Start + u*(yGrid.End-yGrid.Start);
    double z = xGrid.Start + u*(zGrid.End-zGrid.Start);
    fprintf (f1, "%1.12f ", u);
    fprintf (f2, "%1.12f ", u);
    MultiSpline(x,y,z,vals);
    for (int n=0; n<numSplines; n++) {
      fprintf (f1, "%1.16e ", MySplines(n)(x,y,z));
      fprintf (f2, "%1.16e ", vals(n));
    }
    fprintf (f1, "\n");
    fprintf (f2, "\n");
  }
  fclose(f1);
  fclose(f2);
}


void SpeedTest()
{
  int N = 30;
  const int numSplines = 16;
  Array<MyTricubicSpline,1> MySplines(numSplines);
  MultiTricubicSpline MultiSpline;
  LinearGrid xGrid, yGrid, zGrid;
  
  xGrid.Init(0.0, 4.0, N);
  yGrid.Init(0.0, 5.0, N);
  zGrid.Init(0.0, 6.0, N);
  
  Array<double,4> initData(N,N,N,numSplines);
  for (int ix=0; ix<N; ix++)
    for (int iy=0; iy<N; iy++)
      for (int iz=0; iz<N; iz++) 
	for (int n=0; n<numSplines; n++) 
	  initData(ix,iy,iz,n) = drand48();

  for (int n=0; n<numSplines; n++) {
    Array<double,3> temp;
    temp.reference (initData(Range::all(),Range::all(),Range::all(),n));
    MySplines(n).Init (&xGrid, &yGrid, &zGrid, temp);
  }
  MultiSpline.Init (&xGrid, &yGrid, &zGrid, initData);

  Array<double,1> vals(numSplines);
  int numEvals = 16*100000;
  clock_t start, end;

  for (int j=0; j<numSplines; j++)
    vals(j) = MySplines(j)(0.1,0.2,0.3);

  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    for (int j=0; j<numSplines; j++)
      vals(j) = MySplines(j)(x,y,z);
  }
  end = clock();
  fprintf (stderr, "MySplines time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);

  MultiSpline(0.1, 0.2, 0.3, vals);
  start = clock();
  for (int i=0; i<numEvals; i++) {
    double x = xGrid.Start+drand48()*(xGrid.End-xGrid.Start);
    double y = yGrid.Start+drand48()*(yGrid.End-xGrid.Start);
    double z = zGrid.Start+drand48()*(zGrid.End-xGrid.Start);
    MultiSpline(x,y,z,vals);
  }
  end = clock();
  fprintf (stderr, "MultiSpline time = %1.3f sec\n", 
	   (double)(end-start)/CLOCKS_PER_SEC);
}


main()
{
  ValTest();
  SpeedTest();
}
  
