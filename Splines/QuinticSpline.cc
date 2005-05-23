#include "QuinticSpline.h"

#ifdef NOUNDERSCORE 
#define FORT(name) name
#else
#define FORT(name) name ## _
#endif 

extern "C" void FORT(quinat) (int *N, double X[], double Y[], 
			      double B[], double C[], double D[], 
			      double E[], double F[]);




void QuinticSpline::Update()
{
  /// First, use double and triple knots to specify first and second
  /// derivatives at the boundary if we so desire.
  offset=0;
  FX(0) = (*grid)(0);
  FY(0) = Y(0);
  if (!isnan(StartDeriv) && !isinf(StartDeriv)) {
    offset++;
    FX(1) = (*grid)(0);
    FY(1) = StartDeriv;
    if (!isnan(StartDeriv2) && !isinf(StartDeriv2)) { 
      offset++;
      FX(2) = (*grid)(0);
      FY(2) = StartDeriv2;
    }
  }
  
  int i = grid->NumPoints + offset;
  if (!isnan(EndDeriv) && !isinf(EndDeriv)) {
    FX(i) = (*grid)(grid->NumPoints-1);
    FY(i) = EndDeriv;
    i++;
    if (!isnan(EndDeriv2) && !isinf(EndDeriv)) {
      FX(i) = (*grid)(grid->NumPoints-1);
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
  FORT(quinat) (&Fpoints, FX.data(), FY.data(), FB.data(), FC.data(),
		FD.data(), FE.data(), FF.data());

  // Now copy the data into our coefficents
  B(0) = FB(offset);
  C(0) = FC(offset);
  D(0) = FD(offset);
  E(0) = FE(offset);
  F(0) = FF(offset);
  for (int i=1; i<grid->NumPoints; i++) {
    B(i) = FB(i+offset);
    C(i) = FC(i+offset);
    D(i) = FD(i+offset);
    E(i) = FE(i+offset);
  }
  /// Last FF(N-1) is not set by quinat, so we don't copy to avoid a
  /// valgrind error.
  for (int i=1; i<(grid->NumPoints-1); i++)
    F(i) = FF(i+offset);
  F(grid->NumPoints-1) = 0.0;
  UpToDate=true;
}


void QuinticSpline::Write(IOSectionClass &outSection)
{
  outSection.WriteVar("StartDeriv", StartDeriv);
  outSection.WriteVar("EndDeriv", EndDeriv);
  outSection.WriteVar("StartDeriv2", StartDeriv2);
  outSection.WriteVar("EndDeriv2", EndDeriv2);
  outSection.WriteVar("Y", Y);
  
  outSection.NewSection("Grid");
  grid->Write(outSection);
  outSection.CloseSection();
}

void QuinticSpline::Read(IOSectionClass &inSection)
{
  assert(inSection.ReadVar("StartDeriv", StartDeriv));
  assert(inSection.ReadVar("EndDeriv", EndDeriv));
  assert(inSection.ReadVar("StartDeriv2", StartDeriv2));
  assert(inSection.ReadVar("EndDeriv2", EndDeriv2));
  Array<double,1> newY;
  assert(inSection.ReadVar("Y", newY));
  assert(inSection.OpenSection("Grid"));
  Grid *newGrid = ReadGrid(inSection);
  inSection.CloseSection();
  Init (newGrid, newY, StartDeriv, EndDeriv, StartDeriv2, EndDeriv2);
  Update();
}
