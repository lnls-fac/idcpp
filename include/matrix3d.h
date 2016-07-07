#ifndef Matrix3D_H
#define Matrix3D_H

#include <cmath>
#include <ostream>
#include "vector3d.h"

template <typename T = double>
class Matrix3D {
public:
  Matrix3D() : data{0,0,0} {}
  Matrix3D(const Matrix3D<T>& m) : data{m.data[0],m.data[1],m.data[2]} {}
  Matrix3D(const Vector3D<T>& vx, const Vector3D<T>& vy, const Vector3D<T>& vz) : data{vx,vy,vz} {}
  Matrix3D<T>& set_constant(const T& v) { data[0] = Vector3D<T>(v,v,v); data[1] = Vector3D<T>(v,v,v); data[2] = Vector3D<T>(v,v,v); }
  Matrix3D<T>& set_diagonal(const T& v) { data[0].x = v; data[1].y = v; data[2].z = v; return *this; }
  Matrix3D<T>& set_diagonal(const Vector3D<T>& v) { data[0].x = v.x; data[1].y = v.y; data[2].z = v.z; return *this; }
  Vector3D<T>  get_diagonal() const { return Vector3D<T>(data[0].x,data[1].y,data[2].z); }
  Matrix3D<T>& set_rotation_x(const double& angle);
  Matrix3D<T>& set_rotation_y(const double& angle);
  Matrix3D<T>& set_rotation_z(const double& angle);
  Vector3D<T>  row(int i) const { return data[i]; }
  Vector3D<T>  column(int i) const { return Vector3D<T>(data[0][i],data[1][i],data[2][i]); }
  Vector3D<T>  operator *(const Vector3D<T>& v) const;
  Matrix3D<T>  operator *(const Matrix3D<T>& v);
  Matrix3D<T>  operator *(const T& v) { return Matrix3D(v*data[0],v*data[1],v*data[2]); }
  static const Matrix3D<T> I;
  static const Matrix3D<T> rotx90p;
  static const Matrix3D<T> rotx90n;
  static const Matrix3D<T> roty90p;
  static const Matrix3D<T> roty90n;
  static const Matrix3D<T> rotz90p;
  static const Matrix3D<T> rotz90n;
private:
  Vector3D<T> data[3];

};

template <typename T> const Matrix3D<T> Matrix3D<T>::I(Vector3D<T>(1,0,0),Vector3D<T>(0,1,0),Vector3D<T>(0,0,1));
template <typename T> const Matrix3D<T> Matrix3D<T>::rotx90p(Vector3D<T>(+1,0,0),Vector3D<T>(0,0,-1),Vector3D<T>(0,+1,0));
template <typename T> const Matrix3D<T> Matrix3D<T>::rotx90n(Vector3D<T>(+1,0,0),Vector3D<T>(0,0,+1),Vector3D<T>(0,-1,0));
template <typename T> const Matrix3D<T> Matrix3D<T>::roty90p(Vector3D<T>(0,0,+1),Vector3D<T>(0,+1,0),Vector3D<T>(-1,0,0));
template <typename T> const Matrix3D<T> Matrix3D<T>::roty90n(Vector3D<T>(0,0,-1),Vector3D<T>(0,+1,0),Vector3D<T>(+1,0,0));
template <typename T> const Matrix3D<T> Matrix3D<T>::rotz90p(Vector3D<T>(0,-1,0),Vector3D<T>(+1,0,0),Vector3D<T>(0,0,+1));
template <typename T> const Matrix3D<T> Matrix3D<T>::rotz90n(Vector3D<T>(0,+1,0),Vector3D<T>(-1,0,0),Vector3D<T>(0,0,+1));

template <typename T>
Vector3D<T> Matrix3D<T>::operator*(const Vector3D<T>& v) const {
  return Vector3D<T>(data[0].dot(v),data[1].dot(v),data[2].dot(v));
}

template <typename T>
Matrix3D<T> Matrix3D<T>::operator*(const Matrix3D<T>& v) {
  Vector3D<T>  c0 = v.column[0];
  Vector3D<T>  c1 = v.column[1];
  Vector3D<T>  c2 = v.column[2];
  Vector3D<T>& l0 = data[0];
  Vector3D<T>& l1 = data[1];
  Vector3D<T>& l2 = data[2];
  return Matrix3D<T>(Vector3D<T>(l0.dot(c0), l0.dot(c1), l0.dot(c2)),
                    Vector3D<T>(l1.dot(c0), l1.dot(c1), l1.dot(c2)),
                    Vector3D<T>(l2.dot(c0), l2.dot(c1), l2.dot(c2)));
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_x(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  data[0] = Vector3D<T>( 1, 0, 0);
  data[1] = Vector3D<T>( 0,+c,-s);
  data[2] = Vector3D<T>( 0,+s,+c);
  return *this;
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_y(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  data[0] = Vector3D<T>(+c, 0, +s);
  data[1] = Vector3D<T>( 0, 1,  0);
  data[2] = Vector3D<T>(-s, 0, +c);
  return *this;
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_z(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  data[0] = Vector3D<T>(+c, -s, 0);
  data[1] = Vector3D<T>(+s, +c, 0);
  data[2] = Vector3D<T>( 0,  0, 1);
  return *this;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const Matrix3D<T>& v) {
  out << "(" << v.row(0) << "," << v.row(1) << "," << v.row(2) << ")";
  return out;
}


#endif
