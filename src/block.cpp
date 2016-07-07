#include <cmath>
#include <vector>
#include <idmodel.h>

template <typename T>
Matrix3D<T> Block::get_gmatrix(const Vector3D<T>& r) const {
  Matrix3D<T> m;
  for(std::vector<Subblock>::size_type i = 0; i != subblocks.size(); i++) {
    //m = subblocks[i].str * subblocks[i].get_gmatrix(r);
    m = subblocks[i].get_gmatrix(r);
  }
  return m;
}

template <typename T>
Vector3D<T> Block::get_field(const Vector3D<T>& r) const {
  return this->get_gmatrix<T>(r) * this->mag;
}


template <typename T>
Matrix3D<T> Subblock::get_gmatrix(const Vector3D<T>& r) const {
  T x[] = {pos[0] - r[0] - dim[0]/2, pos[0] - r[0] + dim[0]/2};
  T y[] = {pos[1] - r[1] - dim[1]/2, pos[1] - r[1] + dim[1]/2};
  T z[] = {pos[2] - r[2] - dim[2]/2, pos[2] - r[2] + dim[2]/2};
  T gxx,gxy,gxz,gyy,gyz,gzz;
  gxx=gxy=gxz=gyy=gyz=gzz=0;
  for(unsigned int i=0; i<2; i++) {
    for(unsigned int j=0; j<2; j++) {
      for(unsigned int k=0; k<2; k++) {
        int sign = (i+j+k)&1?-1:1;
        T   mod  = std::sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]);
        gxx += sign * std::atan (y[j]*z[k]/(x[i]*mod));
        gyy += sign * std::atan (x[i]*z[k]/(y[j]*mod));
        gzz += sign * std::atan (y[j]*x[i]/(z[k]*mod));
        gxy += - sign * std::log(z[k] + mod);
        gxz += - sign * std::log(y[j] + mod);
        gyz += - sign * std::log(x[i] + mod);
      }
    }
  }
  const double c = 0.25 / M_PI;
  return Matrix3D<T>(Vector3D<T>(c*gxx,c*gxy,c*gxz),
                    Vector3D<T>(c*gxy,c*gyy,c*gyz),
                    Vector3D<T>(c*gxz,c*gyz,c*gzz));
}
