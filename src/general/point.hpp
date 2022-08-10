#ifndef POINT_HPP
#define POINT_HPP

template <typename T>
class Point3 {
public:
  T x;
  T y;
  T z;

  Point3() {};

  Point3(T X, T Y, T Z): x(X), y(Y), z(Z) {};

  ~Point3() {};
};

typedef Point3<double> Point3d;

template <typename T>
inline Point3<T> operator+(const Point3<T>& a, const Point3<T>& b) {
  return Point3<T>(a.x+b.x, a.y+b.y, a.z+b.y);
}

template <typename T, typename S>
inline Point3<T> operator*(const Point3<T>& a, const S& b) {
  return Point3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T, typename S>
inline Point3<T> operator*(const S& b, const Point3<T>& a) {
  return Point3<T>(a.x*b, a.y*b, a.z*b);
}

template <typename T, typename S>
inline Point3<T> operator/(const Point3<T>& a, const S& b) {
  return Point3<T>(a.x/b, a.y/b, a.z/b);
}

#endif
