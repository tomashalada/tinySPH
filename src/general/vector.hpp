#ifndef VECTOR3_HPP
#define VECTOR3_HPP

#include <cmath>
#include "point.hpp"

template <typename T>
class Matrix3;

template <typename T>
class Vector3;

template <typename T>
class Vector3 {
public:
  T x;
  T y;
  T z;

  Vector3() {};
  Vector3(T X, T Y, T Z): x(X), y(Y), z(Z) {};
  Vector3(const Point3<T>& a, const Point3<T>& b): x(b.x-a.x), y(b.y-a.y), z(b.z-a.z) {};

  ~Vector3() {};

  T length() const {
    return sqrt(x*x + y*y + z*z);
  }
};

typedef Vector3<double> Vector3d;

template <typename T>
inline Point3<T> operator+(const Point3<T>& a, const Vector3<T>& b) {
  return Point3<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template <typename T>
inline Point3<T> operator-(const Point3<T>& a, const Vector3<T>& b) {
  return Point3<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template <typename T>
inline Vector3<T> operator-(const Point3<T>& a, const Point3<T>& b) {
  return Vector3<T>(b, a);
}

template <typename T>
inline Vector3<T> operator+(const Vector3<T>&a, const Vector3<T>& b) {
  return Vector3<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template <typename T>
inline Vector3<T> operator-(const Vector3<T>&a, const Vector3<T>& b) {
  return Vector3<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template <typename T, typename S>
inline Vector3<T> operator*(const Vector3<T>& a, const S& b) {
  return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

/*
template <typename T>
inline Vector3<T> operator*(const Vector3<T>& a, const int& b) {
  return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T>
inline Vector3<T> operator*(const Vector3<T>& a, const float& b) {
  return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T>
inline Vector3<T> operator*(const Vector3<T>& a, const double& b) {
  return Vector3<T>(a.x*b, a.y*b, a.z*b);
}
*/


template <typename T, typename S>
inline Vector3<T> operator*(const S& b, const Vector3<T>& a) {
  return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

/*
template <typename T >
inline Vector3<T> operator*(const int& b, const Vector3<T>& a) {
    return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T >
inline Vector3<T> operator*(const float& b, const Vector3<T>& a) {
    return Vector3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T >
inline Vector3<T> operator*(const double& b, const Vector3<T>& a) {
    return Vector3<T>(a.x*b, a.y*b, a.z*b);
}
*/


template <typename T>
inline Vector3<T> operator*(const Vector3<T>& a, const Vector3<T>& b) {
  return Vector3<T>(a.x*b.x, a.y*b.y, a.z*b.z);
}

template <typename T, typename S>
inline Vector3<T> operator/(const Vector3<T>& a, const S& b) {
  return Vector3<T>(a.x/b, a.y/b, a.z/b);
}

template <typename T>
inline Vector3<T> operator+=(Vector3<T>&a, const Vector3<T>& b) {
  a.x+=b.x; a.y+=b.y; a.z+=b.z;
  return a;
}

template <typename T>
inline Vector3<T> operator-=(Vector3<T>&a, const Vector3<T>& b) {
  a.x-=b.x; a.y-=b.y; a.z-=b.z;
  return a;
}

template <typename T, typename S>
inline Vector3<T> operator/=(Vector3<T>&a, const S& b) {
  a.x/=b; a.y/=b; a.z/=b;
  return a;
}

template <typename T>
T dot(const Vector3<T>& a, const Vector3<T>& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

template <typename T>
Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b) {
  return Vector3<T>(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

#endif
