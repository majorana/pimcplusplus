#include "Grid.h"

extern "C" void mktricubw_(double x[], int *nx,
			   double y[], int *ny,
			   double z[], int *nz,
			   double *f, int *nf2, int *nf3,
			   int *ibcxmin, double *bcxmin, 
			   int *ibcxmax, double *bcxmax, int *inb1x,
			   int *ibcymin, double *bcymin, 
			   int *ibcymax, double *bcymax, int *inb1y,
			   int *ibczmin, double *bczmin,
			   int *ibczmax, double *bczmax, int *inb1z,
			   double *wk, int *nwk, 
			   int *ilinx, int *iliny, int *ilinz, int *ier);

extern "C" void evtricub_(double *xget, double *yget, double *zget,
			  double *x, int *nx,
			  double *y, int *ny,
			  double *z, int *nz,
			  int *ilinx, int *iliny, int *ilinz,
			  double *f, int *inf2, int *inf3,
			  int *iselect, double *fval,
			  int *ier);

class TricubicSpline
{
private:
  Array<double,4> f;
  Array<double,1> wk;
  bool UpToDate;
  int Nx, Ny, Nz;
  int xIsLin, yIsLin, zIsLin;
  int IxMinBC, IxMaxBC, IyMinBC, IyMaxBC, IzMinBC, IzMaxBC;
  double dummyDouble;
  int dummyInt;
  int Nwk;
public:
  Grid *xGrid, *yGrid, *zGrid;
  inline double operator()(int ix, int iy, int iz) const
  {
    return (f(iz,iy,iz,0));
  }
  inline double& operator()(int ix, int iy, int iz) 
  {
    UpToDate=false;
    return (f(iz,iy,iz,0));
  }
  inline double operator()(double x, double y, double z); 
  inline void Update();

  inline void Init (Grid *xGrid_, Grid *yGrid_, Grid *zGrid_,
		    Array<double,3> &f_);

  TricubicSpline(Grid *xGrid_, Grid *yGrid_, Grid *zGrid_,
		 Array<double,3> &f_)
  {
    Init(xGrid_, yGrid_, zGrid_, f_);
  }
  TricubicSpline() {};

};




inline double TricubicSpline::operator()(double x, double y, double z) 
{
  if (!UpToDate)
    Update();

  double fval[10];
  int iselect[10];
  iselect[0]=1;
  for (int i=1; i<10; i++)
    iselect[i] = 0;
  int errorCode;

  evtricub_(&x, &y, &z, 
	    xGrid->Points(), &Nx, 
	    yGrid->Points(), &Ny, 
	    zGrid->Points(), &Nz,
	    &xIsLin, &yIsLin, &zIsLin,
	    f.data(), &Nx, &Ny, 
	    iselect, fval, &errorCode);
  if (errorCode != 0 || isnan(fval[0])){
    cerr << "x = " << x << endl;
    cerr << "y = " << y << endl;
    cerr << "z = " << z << endl;
    cerr << "Error in TricubicSpline::operator().\n";
  }

  return (fval[0]);
}




inline void TricubicSpline::Init (Grid *xGrid_, Grid *yGrid_, Grid *zGrid_,
				  Array<double,3> &f_)
{
  xGrid=xGrid_; yGrid=yGrid_; zGrid=zGrid_;
  Nx = xGrid->NumPoints; Ny=yGrid->NumPoints; Nz=zGrid->NumPoints;
  f.resize(Nz, Ny, Nx, 8);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	f(iz, iy, ix, 0) = f_(ix, iy, iz);
  Nwk = 80*Nx*Ny*Nz;
  wk.resize(Nwk);
  IxMinBC=6; IxMaxBC=6;
  IyMinBC=6; IyMaxBC=6;
  IzMinBC=6; IzMaxBC=6;
  xIsLin = (xGrid->Type() == LINEAR);
  yIsLin = (yGrid->Type() == LINEAR);
  zIsLin = (zGrid->Type() == LINEAR);
  Update();
}




inline void TricubicSpline::Update()
{
  int errorCode;
//   cerr << "Nx = " << Nx << endl;
//   cerr << "Ny = " << Ny << endl;
//   cerr << "Nz = " << Nz << endl;
//   for (int i=0; i<Nx; i++)
//     cerr << (*xGrid)(i) << endl;
//   for (int i=0; i<Ny; i++)
//     cerr << (*yGrid)(i) << endl;
//   for (int i=0; i<Nz; i++)
//     cerr << (*zGrid)(i) << endl;


  mktricubw_(xGrid->Points(), &Nx,
	     yGrid->Points(), &Ny,
	     zGrid->Points(), &Nz,
	     f.data(), &Nx, &Ny,
	     &IxMinBC, &dummyDouble, &IxMaxBC, &dummyDouble, &dummyInt,
	     &IyMinBC, &dummyDouble, &IyMaxBC, &dummyDouble, &dummyInt,
	     &IzMinBC, &dummyDouble, &IzMaxBC, &dummyDouble, &dummyInt,
	     wk.data(), &Nwk, &xIsLin, &yIsLin, &zIsLin, &errorCode);
  for (int ix=0; ix<Nx; ix++)
    for (int iy=0; iy<Ny; iy++)
      for (int iz=0; iz<Nz; iz++)
	for (int k=0; k<8; k++)
	  if (isnan(f(iz,iy,ix,k)))
	    fprintf (stderr, "Error at (x,y,z,k) = (%d,%d,%d,%d)\n", 
		     ix, iy, iz, k); 
  UpToDate = true;
}
	     
	     

