#ifndef QUINTIC_SPLINE_H
#define QUINTIC_SPLINE_H

#include "Grid.h"
#include "../IO/InputOutput.h"
#include <iostream>


/// The QuinticSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class QuinticSpline
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  int UpToDate;
  /// The function values on the grid points.
  Array<double, 1> Y;
  /// The parameters to the fortran routine QUINAT
  Array<double, 1> FX, FY, FB, FC, FD, FE, FF;
  int offset;
  /// The start and end first and second derivatives of the function.
  /// If these are NAN's, then we assume natural boundary conditions.
  double StartDeriv, StartDeriv2, EndDeriv, EndDeriv2;
  /// The coefficients of the 5th-order polynomial
  Array<double,1>  B, C, D, E, F;
  double I, J;
  inline GetIJ(double x);
public:
  int NumParams;
  Grid *grid;
  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

  /// Returns the interpolated value.
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double Deriv(double x);
  /// Returns the interpolated second derivative.
  inline double Deriv2(double x);
  /// Returns the interpolated third derivative.
  inline double Deriv3(double x);
  /// Returns the interpolated fourth derivative.
  inline double Deriv4(double x);
  /// Recompute the quintic polynomial coefficients
  void Update();
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(Grid *NewGrid, Array<double,1> NewYs,
		   double startderiv, double startderiv2,
		   double endderiv, double endderiv)
  {
    StartDeriv = startderiv; StartDeriv2 = startderiv2;
    EndDeriv = endderiv;     EndDeriv2 = endderiv2;
    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in QuinticSpline.\n";
	exit(1);
      }
    grid = NewGrid;
    NumParams = grid->NumPoints;
    y.resize(grid->NumPoints);
    d2y.resize(grid->NumPoints);
    y = NewYs;
    Update();
  }

  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init (Grid *NewGrid, Array<double,1> NewYs)
  {
    Init (NewGrid, NewYs, NAN, NAN, NAN, NAN);
  }
  
  /// Simplified constructor.
  inline QuinticSpline (Grid *NewGrid, Array<double,1> NewYs)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init (NewGrid, NewYs, NAN, NAN, NAN, NAN);
  }

  /// Full constructor.
  inline QuinticSpline (Grid *NewGrid, Array<double,1> NewYs,
		      double startderiv, double startderiv2,
		      double endderiv, double endderiv2)
  {
    Init (NewGrid, NewYs, startderiv, startdervi2, 
	  endderiv, endderiv2);
    Update();
  }

  /// Returns the value of the function at the ith grid point.
  inline double Params(int i) const
  {
    return (y(i));
  }
  /// Returns a reference to the value at the ith grid point.
  inline double & Params(int i)
  {
    UpToDate = 0;
    return (y(i));
  }
  void Write(OutputSectionClass &outSection)
  {
    outSection.WriteVar("StartDeriv", StartDeriv);
    outSection.WriteVar("EndDeriv", EndDeriv);
    outSection.WriteVar("y", y);

    outSection.OpenSection("Grid");
    grid->Write(outSection);
    outSection.CloseSection();
  }
  void Read(InputSectionClass &inSection)
  {
    assert(inSection.ReadVar("StartDeriv", StartDeriv));
    assert(inSection.ReadVar("EndDeriv", EndDeriv));
    assert(inSection.ReadVar("y", y));
    NumParams = y.size();
    d2y.resize(NumParams);
    assert(inSection.OpenSection("Grid"));
    grid = ReadGrid(inSection);
    inSection.CloseSection();
    Update();
  }


  /// Trivial constructor
  QuinticSpline()
  {
    UpToDate = 0;
  }
};




inline void GetIJ(double x)
{
  Grid &X = *grid;
#ifdef BZ_DEBUG
  if (x > X.End)
    {
      if (x < (X.End * 1.000000001))
	x = X.End;
      else
	{
	  cerr << "x outside grid in QuinticSpline.\n";
	  cerr << "x = " << x << " X.End = " << X.End << "\n";
	  exit(1);
	}
    }
#endif
  int J = X.ReverseMap(x)+1;
  int I = J-1;
  if (I<0)
    {
      I = 0;
      J = 1;
    }
  if (J>(X.NumPoints-1))
    {
      J = (X.NumPoints-1);
      I = hi-1;
    }
}



inline double QuinticSpline::operator()(double x)
{
  if (!UpToDate)
    Update();

  GetIJ();

  double P = (x-X(I));

  double S = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I);
  return (S);
}



inline double QuinticSpline::Deriv(double x)
{
  if (!UpToDate)
    Update();

  GetIJ();
  double P = x-X(I);

  double S = (((5.0*F(I)*P+4.0*E(I))*P+3.0*D(I))*P+2.0*C(I))*P+B(I);
  return (S);
}



inline double QuinticSpline::Deriv2(double x)
{
  if (!UpToDate)
    Update();

  GetIJ();
  double P = (x-X(I));

  double S = ((20.0*F(I)*P+12.0*E(I))*P+6.0*D(I))*P+2.0*C(I);
  return (S);
}




inline double QuinticSpline::Deriv3(double x)
{
  if (!UpToDate)
    Update();

  GetIJ();
  double P = (x-X(I));
  double S = (60.0*F(I)*P+24.0*E(I))*P+6.0*D(I);
  return (S);
}


inline double QuinticSpline::Deriv4(double x)
{
  if (!UpToDate)
    Update();

  GetIJ();
  double P = (x-X(I));
  double S = 120.0*F(I)*P+24.0*E(I);
  return (S);
}




#endif
