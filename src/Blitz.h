#ifndef MYBLITZ_H
#define MYBLITZ_H

#include <blitz/array.h>
typedef double scalar;

#define NDIM 3

using namespace blitz;
typedef TinyVector<scalar,2> Vec2;
typedef TinyVector<scalar,3> Vec3;
typedef TinyMatrix<scalar,2,2> Mat2;
typedef TinyMatrix<scalar,3,3> Mat3;

typedef TinyVector<scalar,NDIM> dVec;
typedef TinyVector<int,NDIM> dVecInt;

template <class T, int size>
inline TinyVector<T,size> operator-(TinyVector<T,size> v)
{
  TinyVector<T,size> minusv;
  for (int i=0; i<size; i++)
    minusv[i] = -v[i];
  return (minusv);
}

inline Vec2 operator*(const Vec2 &v, scalar s)
{
  Vec2 result;
  result[0] = s*v[0];
  result[1] = s*v[1];
  return (result);
}

inline Vec2 operator*(scalar s, const Vec2 &v)
{
  return (v * s);
}

inline Vec2 operator+(const Vec2 &v1, const Vec2 &v2)
{
  Vec2 result;
  result[0] = v1[0]+v2[0];
  result[1] = v1[1]+v2[1];
  return (result);
}


inline Vec2 operator-(const Vec2 &v1, const Vec2 &v2)
{
  Vec2 result;
  result[0] = v1[0]-v2[0];
  result[1] = v1[1]-v2[1];
  return (result);
}


inline Vec3 operator*(scalar s, const Vec3 &v)
{
  Vec3 result;
  result[0] = s*v[0];
  result[1] = s*v[1];
  result[2] = s*v[2];
  return (result);
}

inline Vec3 operator*(const Vec3 &v, scalar s)
{
  Vec3 result;
  result[0] = s*v[0];
  result[1] = s*v[1];
  result[2] = s*v[2];
  return (result);
}


inline Vec3 operator+(const Vec3 &v1, const Vec3 &v2)
{
  Vec3 result;
  result[0] = v1[0]+v2[0];
  result[1] = v1[1]+v2[1];
  result[2] = v1[2]+v2[2];
  return (result);
}

inline Vec3 operator-(const Vec3 &v1, const Vec3 &v2)
{
  Vec3 result;
  result[0] = v1[0]-v2[0];
  result[1] = v1[1]-v2[1];
  result[2] = v1[2]-v2[2];
  return (result);
}

inline Mat3 operator*(scalar s, const Mat3 &M)
{
  Mat3 sM;
  sM(0,0)=s*M(0,0); sM(0,1)=s*M(0,1); sM(0,2)=s*M(0,2);
  sM(1,0)=s*M(1,0); sM(1,1)=s*M(1,1); sM(1,2)=s*M(1,2);
  sM(2,0)=s*M(2,0); sM(2,1)=s*M(2,1); sM(2,2)=s*M(2,2);
  return (sM);
}

inline double distSqrd(Vec2 a,Vec2 b)
{
  return dot(a-b,a-b);
}

inline double distSqrd(Vec3 a,Vec3 b)
{
  return dot(a-b,a-b);
}



#endif
