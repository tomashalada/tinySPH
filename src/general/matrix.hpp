#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include "point.hpp"
#include "vector.hpp"

template <typename T>
class Matrix3 {
public:
  T a11; T a12; T a13;
  T a21; T a22; T a23;
  T a31; T a32; T a33;

  Matrix3() {};
  Matrix3(T A11, T A12, T A13,
          T A21, T A22, T A23,
          T A31, T A32, T A33)
          : a11(A11), a12(A12), a13(A13),
            a21(A21), a22(A22), a23(A23),
            a31(A31), a32(A32), a33(A33) {};

  ~Matrix3() {};

  T det() const {
    return (a11 * a22 * a33 + \
            a12 * a23 * a31 + \
            a13 * a21 * a32 - \
            a31 * a22 * a13 - \
            a32 * a23 * a11 - \
            a33 * a21 * a12);
  }

};

typedef Matrix3<double> Matrix3d;

/*
template <typename T>
inline Point3<T> operator+(const Point3<T>& a, const Matrix3<T>& b) {
  return Point3<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template <typename T>
inline Point3<T> operator-(const Point3<T>& a, const Matrix3<T>& b) {
  return Point3<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template <typename T>
inline Matrix3<T> operator-(const Point3<T>& a, const Point3<T>& b) {
  return Matrix3<T>(b, a);
}
*/

template <typename T>
inline Matrix3<T> operator+(const Matrix3<T>& a, const Matrix3<T>& b) {
  return Matrix3<T>(a.a11+b.a11, a.a12+b.a12, a.a13+b.a13,
                    a.a21+b.a21, a.a22+b.a22, a.a23+b.a23,
                    a.a31+b.a31, a.a32+b.a32, a.a33+b.a33);
}

template <typename T>
inline Matrix3<T> operator-(const Matrix3<T>& a, const Matrix3<T>& b) {
  return Matrix3<T>(a.a11-b.a11, a.a12-b.a12, a.a13-b.a13,
                    a.a21-b.a21, a.a22-b.a22, a.a23-b.a23,
                    a.a31-b.a31, a.a32-b.a32, a.a33-b.a33);
}

template <typename T, typename S>
inline Matrix3<T> operator*(const Matrix3<T>& a, const S& b) {
  return Matrix3<T>(a.a11*b, a.a12*b, a.a13*b,
                    a.a21*b, a.a22*b, a.a23*b,
                    a.a31*b, a.a32*b, a.a33*b);
}

template <typename T, typename S>
inline Matrix3<T> operator*(const S& b, const Matrix3<T>& a) {
  return Matrix3<T>(a.a11*b, a.a12*b, a.a13*b,
                    a.a21*b, a.a22*b, a.a23*b,
                    a.a31*b, a.a32*b, a.a33*b);
}

template <typename T>
inline Vector3<T> operator*(const Matrix3<T>& a, const Vector3<T>& b) {
  return Vector3<T>(a.a11*b.x + a.a12*b.y + a.a13*b.z,
                    a.a21*b.x + a.a22*b.y + a.a23*b.z,
                    a.a31*b.x + a.a32*b.y + a.a33*b.z);
}

template <typename T>
inline Matrix3<T> operator*(const Matrix3<T>& a, const Matrix3<T>& b) {
  return Matrix3<T>(a.a11*b.a11 + a.a12*b.a21 + a.a13*b.a31,    //c11
                    a.a11*b.a12 + a.a12*b.a22 + a.a13*b.a32,    //c12
                    a.a11*b.a13 + a.a12*b.a23 + a.a13*b.a33,    //c13
                    a.a21*b.a11 + a.a22*b.a21 + a.a23*b.a31,    //c21
                    a.a21*b.a12 + a.a22*b.a22 + a.a23*b.a32,    //c22
                    a.a21*b.a13 + a.a22*b.a23 + a.a23*b.a33,    //c23
                    a.a31*b.a11 + a.a32*b.a21 + a.a33*b.a31,    //c31
                    a.a31*b.a12 + a.a32*b.a22 + a.a33*b.a32,    //c32
                    a.a31*b.a13 + a.a32*b.a23 + a.a33*b.a33);   //c33

}

template <typename T, typename S>
inline Matrix3<T> operator/(const Matrix3<T>& a, const S& b) {
  return Matrix3<T>(a.a11/b, a.a12/b, a.a13/b,
                    a.a21/b, a.a22/b, a.a23/b,
                    a.a31/b, a.a32/b, a.a33/b);
}

template <typename T>
inline Matrix3<T> operator+=(Matrix3<T>&a, const Matrix3<T>& b) {
  a.a11+=b.a11; a.a12+=b.a12; a.a13+=b.a13;
  a.a21+=b.a21; a.a22+=b.a22; a.a23+=b.a23;
  a.a31+=b.a31; a.a32+=b.a32; a.a33+=b.a33;
  return a;
}

template <typename T>
inline Matrix3<T> operator-=(Matrix3<T>&a, const Matrix3<T>& b) {
  a.a11-=b.a11; a.a12-=b.a12; a.a13-=b.a13;
  a.a21-=b.a21; a.a22-=b.a22; a.a23-=b.a23;
  a.a31-=b.a31; a.a32-=b.a32; a.a33-=b.a33;
  return a;
}

template <typename T, typename S>
inline Matrix3<T> operator/=(Matrix3<T>&a, const S& b) {
  a.a11/=b.a11; a.a12/=b.a12; a.a13/=b.a13;
  a.a21/=b.a21; a.a22/=b.a22; a.a23/=b.a23;
  a.a31/=b.a31; a.a32/=b.a32; a.a33/=b.a33;
  return a;
}

template <typename T>
Matrix3<T> inv(const Matrix3<T>& a) {
  double det = a.det();
  return Matrix3<T>( (a.a22*a.a33-a.a23*a.a32)/det,
                    -(a.a12*a.a33-a.a13*a.a32)/det,
                     (a.a12*a.a23-a.a13*a.a22)/det,
                    -(a.a21*a.a33-a.a23*a.a31)/det,
                     (a.a11*a.a33-a.a13*a.a31)/det,
                    -(a.a11*a.a23-a.a13*a.a21)/det,
                     (a.a21*a.a32-a.a22*a.a31)/det,
                    -(a.a11*a.a32-a.a12*a.a31)/det,
                     (a.a11*a.a22-a.a12*a.a21)/det);
}

template <typename T>
Vector3<T> solve(const Matrix3<T>& a, const Vector3<T>& b) {
  return inv(a)*b;
}

#endif
