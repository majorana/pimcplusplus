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

template <class T, int size>
inline bool operator==(TinyVector<T,size> v1,TinyVector<T,size> v2)
{
  bool equals=true;


  for (int i=0; i<size; i++)
    equals = (v1[i]==v2[i]) && equals;
  return (equals);
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

inline Mat3 operator*(const Mat3 &A, const Mat3 &B)
{
  Mat3 AB;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      AB(i,j) = 0.0;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++)
	AB(i,j) += A(i,k)*B(k,j);
  return (AB);
}

inline Mat3 operator+(const Mat3 &A, const Mat3 &B)
{
  Mat3 ApB;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      ApB(i,j) = A(i,j) + B(i,j);
  return (ApB);
}

inline Vec3 operator*(const Mat3& A, const Vec3 &v)
{
  Vec3 Av;
  Av[0] = A(0,0)*v[0] + A(0,1)*v[1] + A(0,2)*v[2];
  Av[1] = A(1,0)*v[0] + A(1,1)*v[1] + A(1,2)*v[2];
  Av[2] = A(2,0)*v[0] + A(2,1)*v[1] + A(2,2)*v[2];
  return Av;
}

inline double operator*(const Vec3 &v1, const Vec3 &v2)
{
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}


inline double distSqrd(Vec2 a,Vec2 b)
{
  return dot(a-b,a-b);
}

inline double distSqrd(Vec3 a,Vec3 b)
{
  return dot(a-b,a-b);
}


template<class T> class SymmArray
{
private:
  Array<T,1> A;
  int N;
  inline int index(int row, int col) const 
  { return ((row > col) ? ((row*(row+1)>>1)+col) : ((col*(col+1)>>1)+row)); }
public:

  inline void resize(int n) 
  {
    N = (n*(n+1))>>1;
    A.resize(N);
  }

  inline int rows() const
  { return N; }
  
  inline T operator()(int i, int j) const
  { return (A(index(i,j))); }

  inline T& operator()(int i, int j)
  { return (A(index(i,j))); }

  inline SymmArray<T> (const SymmArray<T> &B)
  {
    resize(B.N);
    A=B.A;
  }
  inline SymmArray<T>& operator=(const SymmArray<T>& B)
  {
    A=B.A;
  }
  inline SymmArray<T>()
  { N=0; }
};




#endif
