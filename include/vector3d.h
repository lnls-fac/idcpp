#ifndef Vector3D_H
#define Vector3D_H

#include <iostream>
#include <ostream>
#include <cmath>

template <typename T = double>
class Vector3D {

public:
  Vector3D() : x(0), y(0), z(0) {}
  Vector3D(const T& x_, const T& y_=0, const T& z_=0) : x(x_), y(y_), z(z_) {}
  Vector3D(const Vector3D& v) : x(v.x), y(v.y), z(v.z) {}
  T x,y,z;
  Vector3D<T> operator  -  () const;
  Vector3D<T> operator  +  (const Vector3D<T>& v) const;
  Vector3D<T> operator  -  (const Vector3D<T>& v) const;
  Vector3D<T> operator  *  (const Vector3D<T>& v) const;
  Vector3D<T> operator  *  (const T& t) const;
  Vector3D<T> operator  /  (const T& t) const;
  Vector3D<T>& operator *= (const T& t);
  Vector3D<T>& operator /= (const T& t);
  Vector3D<T>& operator += (const Vector3D<T>& v) { x += v.x; y += v.y; z += v.z; return *this; }
  Vector3D<T>& operator -= (const Vector3D<T>& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
  bool operator != (const Vector3D& v) const { return (x!=v.x)or(y!=v.y)or(z!=v.z); }
  bool operator == (const Vector3D& v) const { return (x==v.x)and(y==v.y)and(z==v.z); }
  inline T dot(const Vector3D<T> v) const { return this->x*v.x + this->y*v.y + this->z*v.z; }

private:
  T data;
};

template <typename T>
Vector3D<T> Vector3D<T>::operator+(const Vector3D<T>& v) const {
  return Vector3D<T>(this->x+v.x,this->y+v.y,this->z+v.z);
}

template <typename T>
Vector3D<T> Vector3D<T>::operator-(const Vector3D<T>& v) const {
  return Vector3D<T>(this->x-v.x,this->y-v.y,this->z-v.z);
}

template <typename T>
Vector3D<T> Vector3D<T>::operator*(const Vector3D<T>& v) const {
  return Vector3D<T>(this->x*v.x,this->y*v.y,this->z*v.z);
}

template <typename T>
Vector3D<T> Vector3D<T>::operator-() const {
  return Vector3D<T>(-this->x,-this->y,-this->z);
}

template <typename T>
Vector3D<T> Vector3D<T>::operator*(const T& t) const {
  return Vector3D<T>(this->x*t,this->y*t,this->z*t);
}

template <typename T>
Vector3D<T> Vector3D<T>::operator/(const T& t) const {
  return Vector3D<T>(this->x/t,this->y/t,this->z/t);
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator*=(const T& t) {
  this->x *= t; this->y *= t; this->z *= t;
  return *this;
}

template <typename T>
Vector3D<T>& Vector3D<T>::operator/=(const T& t) {
  this->x /= t; this->y /= t; this->z /= t;
  return *this;
}

template <typename T>
Vector3D<T> operator* (const int& v1, const Vector3D<T>& v2) {
  return Vector3D<T>(v1*v2.x,v1*v2.y,v1*v2.z);
}

template <typename T>
std::ostream& operator <<(std::ostream& out, const Vector3D<T>& v) {
  out << "(" << v.x << "," << v.y << "," << v.z << ")";
  return out;
}

#endif
