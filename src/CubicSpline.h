#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "Grid.h"

class CubicSpline
{
private:
  int UpToDate;
  Array<double, 1> y;   // The values of the function at the data points
  Array<double, 1> d2y;  // The second derivatives of the function
public:
  int NumParams;
  Grid *grid;
  double StartDeriv, EndDeriv;

  double operator()(double x);
  double Deriv(double x);
  double Deriv2(double x);
  double Deriv3(double x);
  void Update();
  
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

  inline void Init (Grid *NewGrid, Array<double,1> NewYs)
  {
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }

  inline CubicSpline (Grid *NewGrid, Array<double,1> NewYs)
  {
    StartDeriv = EndDeriv = 5.0e30;
    Init (NewGrid, NewYs, 5.0e30, 5.0e30);
  }

  inline CubicSpline (Grid *NewGrid, Array<double,1> NewYs,
		      double startderiv, double endderiv)
  {
    Init (NewGrid, NewYs, startderiv, endderiv);
    Update();
  }

  inline double Params(int i) const
  {
    return (y(i));
  }
  inline double & Params(int i)
  {
    UpToDate = 0;
    return (y(i));
  }

  CubicSpline()
  {
    UpToDate = 0;
  }
};


class MultiCubicSpline
{
private:
  Array<bool,1> UpToDate;
  Array<double, 2> y;   // The values of the function at the data points
  Array<double, 2> d2y;  // The second derivatives of the function
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
