#include "GracePlot.h"

#include <unistd.h>

main()
{
  Array<double,1> x(100), y(100);
  for (int i=0; i<100; i++) {
    x(i) = 4.0*i*M_PI/99.0;
    y(i) = sin(x(i));
  }
  GracePlot plot1; 
  plot1.AddCurve();
  plot1.Curves(0).SetData(x,y);
  plot1.Curves(0).SetColor ("red");
  plot1.Curves(0).SetLineWidth(4.0);
  //  plot1.SetBounds(0.0, -1.0, 4.0*M_PI, 1.0);
  
  sleep(5);
  GracePlot plot2;
  for (int i=0; i<100; i++) {
    x(i) = 4.0*i*M_PI/99.0;
    y(i) = cos(x(i));
  }
  plot2.AddCurve();
  plot2.Curves(0).SetData(x,y);
  plot2.Curves(0).SetColor("blue");
  plot2.Curves(0).SetLineWidth(2.0);
  //  plot2.SetBounds(0.0, -1.0, 4.0*M_PI, 1.0);
  plot2.Redraw();
  sleep(5);

  plot1.MakeCurrent();
  sleep(5);
}
