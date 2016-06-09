#ifndef FIELDMAP_VECTOR3D_H
#define FIELDMAP_VECTOR3D_H

#include <ostream>

template <typename T = double>
class Vector3D {
public:
	Vector3D() : x(0), y(0), z(0) {}
	Vector3D(const T& x_, const T& y_, const T& z_) : x(x_), y(y_), z(z_) {}
	Vector3D(const Vector3D& v) : x(v.x), y(v.y), z(v.z) {}
	T x,y,z;
	Vector3D operator*(const T& v){return Vector3D(this->x*v, this->y*v, this->z*v);}
	Vector3D operator+(const Vector3D & v){return Vector3D(this->x+v.x, this->y+v.y, this->z+v.z);}
	Vector3D operator-(const Vector3D & v){return Vector3D(this->x-v.x, this->y-v.y, this->z-v.z);}
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vector3D<T>& v)
{
  os << "(" << v.x << "," << v.y << "," << v.z << ")";
  return os;
}

template <typename T>
Vector3D<T> operator*(const T& v1, const Vector3D<T> & v2){
	return Vector3D<T> (v2.x*v1, v2.y*v1, v2.z*v1);
}

#endif
