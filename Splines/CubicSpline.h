#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "Grid.h"
#include "../IO/InputOutput.h"
#include <iostream>


/// The CubicSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid. 
class CubicSpline
{
private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  int UpToDate;
  /// The function values on the grid points.
  Array<double, 1> y;   
  /// The second derivatives of the function
  Array<double, 1> d2y;  
public:
  int NumParams;
  Grid *grid;
  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

  /// Returns the interpolated value.
  double operator()(double x);
  /// Returns the interpolated first derivative.
  double Deriv(double x);
  /// Returns the interpolated second derivative.
  double Deriv2(double x);
  /// Returns the interpolated third derivative.
  double Deriv3(double x);
  /// Recompute the second derivatives from the function values
  void Update();
  
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(Grid *NewGrid, Array<double,1> NewYs,
		   double startderiv, double endderiv)
  {
    StartDeriv = startderiv;
    EndDeriv = endderiv;
    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSpline.\n";
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
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }
  
  /// Simplified constructor.
  inline CubicSpline (Grid *NewGrid, Array<double,1> NewYs)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }

  /// Full constructor.
  inline CubicSpline (Grid *NewGrid, Array<double,1> NewYs,
		      double startderiv, double endderiv)
  {
    Init (NewGrid, NewYs, startderiv, endderiv);
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
  CubicSpline()
  {
    UpToDate = 0;
  }
};


/// The MulitCubicSpline class is nearly identical to the CubicSpline
/// class, except that it stores the values of several functions at
/// the same grid point.  This can be used to simultaneously
/// interpolate several functions at the same point without redundant
/// computation.  
class MultiCubicSpline
{
private:
  /// Stores whether each functions second derivatives are in sync
  /// with the stored function values.
  Array<bool,1> UpToDate;
  Array<double, 2> y;   //< The values of the function at the data points
  Array<double, 2> d2y;  //< The second derivatives of the function
public:
  int NumGridPoints, NumSplines;
  Grid *grid;
  Array<double, 1> StartDeriv, EndDeriv;

  double operator()(int, double x);
  void operator()(double x, Array<double,1> &yVec);
  double Deriv(int i, double x);
  void   Deriv (double x, Array<double,1> &deriv);
  double Deriv2(int i, double x);
  void   Deriv2 (double x, Array<double,1> &deriv2);
  double Deriv3(int i, double x);
  void   Deriv3 (double x, Array<double,1> &deriv3);
  void Update(int i);
  
  inline void Init(Grid *NewGrid, const Array<double,2> &NewYs,
		   const Array<double,1> &startderiv, 
		   const Array<double,1> &endderiv)
  {
    grid = NewGrid;
    NumGridPoints = grid->NumPoints;
    NumSplines = NewYs.cols(); 
    UpToDate.resize(NumSplines);
    StartDeriv.resize(NumSplines);
    StartDeriv = startderiv;
    EndDeriv.resize(NumSplines);
    EndDeriv = endderiv;
    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSpline.\n";
	exit(1);
      }
    
    y.resize(NumGridPoints,NumSplines);
    d2y.resize(NumGridPoints,NumSplines);
    y = NewYs;
    for (int i=0; i<NumSplines; i++)
      Update(i);
  }

  inline void Init (Grid *NewGrid, const Array<double,2> &NewYs)
  {
    grid = NewGrid;
    NumGridPoints = grid->NumPoints;
    NumSplines = NewYs.cols();
    UpToDate.resize(NumSplines);
    StartDeriv.resize(NumSplines);
    EndDeriv.resize(NumSplines);
    StartDeriv = 5.0e30;
    EndDeriv = 5.0e30;    
    if (NewGrid->NumPoints != NewYs.rows())
      {
	cerr << "Size mismatch in CubicSpline.\n";
	exit(1);
      }
   
    y.resize(NumGridPoints,NumSplines);
    d2y.resize(NumGridPoints,NumSplines);
    y = NewYs;
    for (int i=0; i<NumSplines; i++)
      Update(i);
  }

  inline MultiCubicSpline (Grid *NewGrid, Array<double,2> NewYs,
			   Array<double,1> startderiv, 
			   Array<double,1> endderiv)
  {
    Init (NewGrid, NewYs, startderiv, endderiv);
  }

  inline MultiCubicSpline (Grid *NewGrid, Array<double,2> NewYs)
  {
    Init (NewGrid, NewYs);
  }

  inline double Params(int i, int j) const
  {
    return (y(i,j));
  }

  inline double & Params(int i, int j)
  {
    UpToDate(j) = false;
    return (y(i,j));
  }

  MultiCubicSpline()
  {
  }
};



#endif
