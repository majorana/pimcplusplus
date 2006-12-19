#ifndef COLOR_MAP_H
#define COLOR_MAP_H

#include <Common/Splines/CubicBspline.h>

typedef enum {BLUE_WHITE_RED, JET, GRAY} ColorMapType;

class ColorMap
{
private:
  CubicBspline Splines[4];
  bool Initialized;
  double Min, Max;
public:
  inline void operator()(double x, TinyVector<double,3> &val);
  inline void operator()(double x, TinyVector<double,4> &val);
  void Init (double min, double max, ColorMapType map=BLUE_WHITE_RED);
  ColorMap() : Initialized(false)
  {

  }
};

inline void
ColorMap::operator() (double x, TinyVector<double,3> &val)
{
  x = max (Min, x);
  x = min (Max, x);
  val[0] = max(0.0, min(Splines[0](x), 1.0));
  val[1] = max(0.0, min(Splines[1](x), 1.0));
  val[2] = max(0.0, min(Splines[2](x), 1.0));
}

inline void
ColorMap::operator() (double x, TinyVector<double,4> &val)
{
  if (Initialized) {
    x = max (Min, x);
    x = min (Max, x);
    val[0] = max(0.0, min(Splines[0](x), 1.0));
    val[1] = max(0.0, min(Splines[1](x), 1.0));
    val[2] = max(0.0, min(Splines[2](x), 1.0));
    val[3] = max(0.0, min(Splines[3](x), 1.0));
  }
  else
    val = 0.0;
}


#endif
