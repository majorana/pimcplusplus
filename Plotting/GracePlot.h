#ifndef GRACE_PLOT_H
#define GRACE_PLOT_H
#include "../Blitz.h"
#include <grace_np.h>
#include <list>


class GraceCurve
{
private:

public:
  int PlotNum, CurveNum;
  void SetLineWidth(double width);
  void SetColor(string color);
  void SetData(const Array<double,1> &x, const Array<double,1> &y);
  void Redraw();
};

class GracePlot
{
private:
  static list<int> Plots;
  int MyPlot, MyNumCurves;
  void SendCommand (stringstream &cmd);
public:
  void MakeCurrent();
  Array<GraceCurve, 1> Curves;
  inline int NumCurves() { return MyNumCurves; }
  void AddCurve();
  void SetTitle(string title);
  void Redraw();
  void SetBounds (double xmin, double ymin, double xmax, double ymax);

  GracePlot();
  ~GracePlot();
};

#endif
