#ifndef Matrix3D_H
#define Matrix3D_H

#include <cmath>
#include <ostream>
#include "vector3d.hpp"

template <typename T = double>
class Matrix3D {

public:

  Matrix3D() : vx(0.0,0.0,0.0), vy(0.0,0.0,0.0), vz(0.0,0.0,0.0) {}
  Matrix3D(const Matrix3D<T>& m) : vx(m.vx), vy(m.vy), vz(m.vz) {}
  Matrix3D(const Vector3D<T>& _vx, const Vector3D<T>& _vy, const Vector3D<T>& _vz) : vx(_vx), vy(_vy), vz(_vz) {}
  Matrix3D<T>& set_constant(const T& _v) { vx = Vector3D<T>(_v,_v,_v); vy = Vector3D<T>(_v,_v,_v); vz = Vector3D<T>(_v,_v,_v); return *this;}
  Matrix3D<T>& set_diagonal(const T& _v) { vx.x = _v; vy.y = _v; vz.z = _v; return *this; }
  Matrix3D<T>& set_diagonal(const Vector3D<T>& _v) { vx.x = _v.x; vy.y = _v.y; vz.z = _v.z; return *this; }
  Vector3D<T>  get_diagonal() const { return Vector3D<T>(vx.x,vy.y,vz.z); }
  Matrix3D<T>& set_rotation_x(const double& angle);
  Matrix3D<T>& set_rotation_y(const double& angle);
  Matrix3D<T>& set_rotation_z(const double& angle);
  Vector3D<T>  column(int i) const;
  Vector3D<T>  row(int i) const;
  Vector3D<T>  operator *(const Vector3D<T>& v) const;
  Matrix3D<T>  operator *(const Matrix3D<T>& v);
  Matrix3D<T>  operator *(const T& v) { return Matrix3D(v*vx,v*vy,v*vz); }
  static const Matrix3D<T> I() { return Matrix3D<T>(Vector3D<double>(1.0,0.0,0.0),Vector3D<double>(0.0,1.0,0.0),Vector3D<double>(0.0,0.0,1.0));}
  static const Matrix3D<T> rotx90p() {return Matrix3D<T>(Vector3D<T>(+1.0,0.0,0.0),Vector3D<T>(0.0,0.0,-1.0),Vector3D<T>(0.0,+1.0,0.0));}
  static const Matrix3D<T> rotx90n() {return Matrix3D<T>(Vector3D<T>(+1.0,0.0,0.0),Vector3D<T>(0.0,0.0,+1.0),Vector3D<T>(0.0,-1.0,0.0));}
  static const Matrix3D<T> roty90p() {return Matrix3D<T>(Vector3D<T>(0.0,0.0,+1.0),Vector3D<T>(0.0,+1.0,0.0),Vector3D<T>(-1.0,0.0,0.0));}
  static const Matrix3D<T> roty90n() {return Matrix3D<T>(Vector3D<T>(0.0,0.0,-1.0),Vector3D<T>(0.0,+1.0,0.0),Vector3D<T>(+1.0,0.0,0.0));}
  static const Matrix3D<T> rotz90p() {return Matrix3D<T>(Vector3D<T>(0.0,-1.0,0.0),Vector3D<T>(+1.0,0.0,0.0),Vector3D<T>(0.0,0.0,+1.0));}
  static const Matrix3D<T> rotz90n() {return Matrix3D<T>(Vector3D<T>(0.0,+1.0,0.0),Vector3D<T>(-1.0,0.0,0.0),Vector3D<T>(0.0,0.0,+1.0));}

private:

  Vector3D<T> vx;
  Vector3D<T> vy;
  Vector3D<T> vz;

};

template <typename T>
Vector3D<T> Matrix3D<T>::row(int i) const{
  switch (i) {
    case 0: return this->vx;
    case 1: return this->vy;
    case 2: return this->vz;
  };
}

template <typename T>
Vector3D<T> Matrix3D<T>::column(int i) const{
  switch (i) {
    case 0: return Vector3D<T>(this->vx.x,this->vy.x,this->vz.x);
    case 1: return Vector3D<T>(this->vx.y,this->vy.y,this->vz.y);
    case 2: return Vector3D<T>(this->vx.z,this->vy.z,this->vz.z);
  };
}

template <typename T>
Vector3D<T> Matrix3D<T>::operator*(const Vector3D<T>& v) const {
  return Vector3D<T>(vx.dot(v),vy.dot(v),vz.dot(v));
}

template <typename T>
Matrix3D<T> Matrix3D<T>::operator*(const Matrix3D<T>& v) {
  Vector3D<T>  c0 = v.column(0);
  Vector3D<T>  c1 = v.column(1);
  Vector3D<T>  c2 = v.column(2);
  Vector3D<T>& l0 = vx;
  Vector3D<T>& l1 = vy;
  Vector3D<T>& l2 = vz;
  return Matrix3D<T>(Vector3D<T>(l0.dot(c0), l0.dot(c1), l0.dot(c2)),
                    Vector3D<T>(l1.dot(c0), l1.dot(c1), l1.dot(c2)),
                    Vector3D<T>(l2.dot(c0), l2.dot(c1), l2.dot(c2)));
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_x(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  vx = Vector3D<T>( 1.0, 0.0, 0.0);
  vy = Vector3D<T>( 0.0, +c, -s);
  vz = Vector3D<T>( 0.0, +s, +c);
  return *this;
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_y(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  vx = Vector3D<T>(+c, 0.0, +s);
  vy = Vector3D<T>( 0.0, 1.0,  0.0);
  vz = Vector3D<T>(-s, 0.0, +c);
  return *this;
}

template <typename T>
Matrix3D<T>& Matrix3D<T>::set_rotation_z(const double& angle) {
  double c = std::cos(angle);
  double s = std::sin(angle);
  vx = Vector3D<T>(+c, -s, 0.0);
  vy = Vector3D<T>(+s, +c, 0.0);
  vz = Vector3D<T>( 0.0,  0.0, 1.0);
  return *this;
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const Matrix3D<T>& v) {
  out << "(" << v.row(0) << "," << v.row(1) << "," << v.row(2) << ")";
  return out;
}


#endif
