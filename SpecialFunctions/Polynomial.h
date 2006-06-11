/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef POLYNOMIAL_CLASS_H
#define POLYNOMIAL_CLASS_H
#include "../Blitz.h"

/// Simple polynomial class.  []'s refer to coefficients.  
/// ()'s refer to polynomial evaluation.
class PolynomialClass
{
  Array<double,1> C;
public:
  inline double operator[](int i) const
  { return (C(i)); }
  
  inline double& operator[](int i)
  { return (C(i)); }
  
  inline int Order() const
  { return (C.size()-1); }
  
  inline void SetOrder(int order)
  { C.resize(order+1); }
  
  inline double operator()(double x) const
  {
    int order = Order();
    double sum=0.0;
    double x2n=1.0;
    for (int n=0; n<=order; n++) {
      sum += C(n)*x2n;
      x2n*=x;
    }
    return (sum);
  }
  
  inline PolynomialClass& operator=(const PolynomialClass &P)
  {
    C.resize(P.Order()+1);
    for (int i=0; i<=P.Order(); i++)
      C(i) = P[i];
    return (*this);
  }
  inline PolynomialClass& operator+=(const PolynomialClass &P)
  {
    assert(Order() == P.Order());
    for (int i=0; i<=P.Order(); i++)
      C(i) += P[i];
    return (*this);
  }

  inline PolynomialClass& operator=(double x)
  {
    for (int i=0; i<=Order(); i++)
      C(i) = x;
  }

  PolynomialClass (const PolynomialClass &P) 
  {
    C.resize(P.Order()+1);
    for (int i=0; i<=P.Order(); i++)
      C(i) = P[i];
  }

  PolynomialClass (int order)
  {
    C.resize(order+1);
    C = 0.0;
  }
  PolynomialClass() {}
};




inline PolynomialClass operator+(const PolynomialClass &a, 
				 const PolynomialClass &b)
{
  int order = (a.Order()>b.Order()) ? a.Order() : b.Order();
  PolynomialClass c(order);
  for (int i=0; i<=order; i++)
    c[i] = 0.0;
  for (int i=0; i<=a.Order(); i++)
    c[i] += a[i];
  for (int i=0; i<=b.Order(); i++)
    c[i] += b[i];
  return (c);
}

inline PolynomialClass operator-(const PolynomialClass &a, 
				 const PolynomialClass &b)
{
  int order = max(a.Order(),b.Order());
  PolynomialClass c(order);
  for (int i=0; i<=order; i++)
    c[i] = 0.0;
  for (int i=0; i<=a.Order(); i++)
    c[i] += a[i];
  for (int i=0; i<=b.Order(); i++)
    c[i] -= b[i];
  return (c);
}


inline PolynomialClass operator*(const PolynomialClass &a, 
				 const PolynomialClass &b)
{
  int order = a.Order()+b.Order();
  PolynomialClass c(order);
  for (int i=0; i<order; i++)
    c[i] = 0.0;
  for (int ai=0; ai<=a.Order(); ai++) 
    for (int bi=0; bi<=b.Order(); bi++) 
      c[ai+bi] += a[ai]*b[bi];
  return (c);
}


inline PolynomialClass operator*(double a, PolynomialClass &P)
{
  int order = P.Order();
  PolynomialClass aP(order);
  for (int i=0; i<=order; i++)
    aP[i] = a*P[i];
  return (aP);
}

inline PolynomialClass operator*(PolynomialClass &P, double a)
{ return (a*P); }

#endif
