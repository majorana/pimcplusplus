#include "QuinticSpline.h"

extern "C" void quinat_ (int *N, double X[], double Y[], 
			 double B[], double C[], double D[], 
			 double E[] double F[]);



void QuinticSpline::Init(Grid *NewGrid, Array<double,1> NewY,
			 double startderiv, double endderiv,
			 double startderiv2, double enderiv2)
{
  StartDeriv = startderiv; StartDeriv2 = startderiv2;
  EndDeriv = endderiv; EndDeriv2 = endderiv2;
  int NumParams = grid->NumPoints;
  if (!isnan(StartDeriv))
    {
      NumParams++;
      if (!isnan(StartDeriv2))
	NumParams++;
    }
  if (!isnan(EndDeriv))
    {
      NumParams++;
      if (!isnan(EndDeriv2))
	NumParams++;
    }
  FX.resize(NumParams);
  FY.resize(NumParams);
  FB.resize(NumParams);
  FC.resize(NumParams);
  FD.resize(NumParams);
  FE.resize(NumParams);
  FF.resize(NumParams);

  Update();
}

void QuinticSpline::Update()
{

  /// First, use double and triple knots to specify first and second
  /// derivatives at the boundary if we so desire.
  offset=0;
  FX(0) = (*grid)(0);
  FY(0) = NewY(0);
  if (!isnan(StartDeriv)) {
    offset++;
    FX(1) = (*grid)(0);
    FY(1) = StartDeriv;
    if (!isnan(EndDeriv)) {
      offset++;
      FX(2) = (*grid)(0);
      FY(2) = StartDeriv2;
    }
  }
  
  int i = grid->NumPoints + offset;
  if (!(isnan(EndDeriv))) {
    FX(i) = (*grid)(grid->NumPoints-1);
    FY(i) = EndDeriv;
    i++;
    if (!(isnan(EndDeriv2))) {
      FX(i) = *(grid)(grid->NumPoints-1);
      FY(i) = EndDeriv2;
    }
  }

  // Now fill in the rest of the values.
  for (int i=1; i<grid->NumPoints; i++) {
    FX(i+offset) = (*grid)(i);
    FY(i+offset) = Y(i);
  }
  int Fpoints = FX.size();
  // Call FORTRAN routine
  quinat_ (&Fpoints, FX.data(), FY.data(), FB.data(), FC.data(),
	   FD.data(), FE.data(), FF.data());

  // Now copy the data into our coefficents
  Y(0) = FY(0);
  

}
