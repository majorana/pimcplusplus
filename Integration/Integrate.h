#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "../Blitz.h"
#include "../Splines/Grid.h"

void RungeKutta4 (const Grid &grid, 
		  int StartPoint, int EndPoint,
		  Array<scalar, 2> &Result,
		  Array<scalar,1> (*DerivFunc)(scalar x,
					       Array<scalar,1> y,
					       void *Params),
		  void *Params);



void IntegrateFirstOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<scalar,1> &Result,
			   scalar (*DerivFunc)(scalar x, scalar y,
					     void *Params),
			  void *Params);
void IntegrateFirstOrderNS (const Grid &grid,
			    int StartPoint, int EndPoint,
			    Array<scalar,1> &Result,
			    scalar (*DerivFunc)(scalar x, scalar y,
						void *Params),
			    void *Params);
void IntegrateSecondOrder (const Grid &grid,
			   int StartPoint, int EndPoint,
			   Array<Vec2,1> &Result,
			   Vec2 (*DerivFunc)(scalar x, Vec2 y,
					     void *Params),
			   void *Params);


#endif
