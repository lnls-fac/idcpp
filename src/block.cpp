#include <cmath>
#include <vector>
#include <api.h>

SubVolume::SubVolume(const SubVolume &obj){
  this->dim = obj.dim;
  this->pos = obj.pos;
  this->str = obj.str;
}

Matrix3D<double> SubVolume::get_gmatrix(const Vector3D<double>& r) const {
  double x[] = {pos.x - r.x - dim.x/2, pos.x - r.x + dim.x/2};
  double y[] = {pos.y - r.y - dim.y/2, pos.y - r.y + dim.y/2};
  double z[] = {pos.z - r.z - dim.z/2, pos.z - r.z + dim.z/2};
  double gxx,gxy,gxz,gyy,gyz,gzz;
  gxx=gxy=gxz=gyy=gyz=gzz=0;
  for(unsigned int i=0; i<2; i++) {
    for(unsigned int j=0; j<2; j++) {
      for(unsigned int k=0; k<2; k++) {
        int sign = (i+j+k)&1?-1:1;
        double   mod  = std::sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]);
        gxx += sign * std::atan (y[j]*z[k]/(x[i]*mod));
        gyy += sign * std::atan (x[i]*z[k]/(y[j]*mod));
        gzz += sign * std::atan (y[j]*x[i]/(z[k]*mod));
        gxy += - sign * std::log(z[k] + mod);
        gxz += - sign * std::log(y[j] + mod);
        gyz += - sign * std::log(x[i] + mod);
      }
    }
  }
   double c = 0.25 / M_PI;
  return Matrix3D<double>(Vector3D<double>(c*gxx,c*gxy,c*gxz),
                    Vector3D<double>(c*gxy,c*gyy,c*gyz),
                    Vector3D<double>(c*gxz,c*gyz,c*gzz));
}


Block::Block(const Block &obj){
  this->mag = obj.mag;
  this->subvolumes = obj.subvolumes;
}

Matrix3D<double> Block::get_gmatrix(const Vector3D<double>& r) const {
  Matrix3D<double> m;
  for(std::vector<SubVolume>::size_type i = 0; i != subvolumes.size(); i++) {
    m = subvolumes[i].get_gmatrix(r);
  }
  return m;
}


Vector3D<double> Block::field(const Vector3D<double>& r) const {
  return this->get_gmatrix(r) * this->mag;
}

void Block::set_mag(Vector3D<double> mag_){
  this->mag = mag_;
}

void Block::set_pos(Vector3D<double> pos_){
  this->subvolumes[0].pos = pos_;
}

void Block::set_dim(Vector3D<double> dim_){
  this->subvolumes[0].dim = dim_;
}

void Block::set_subvolume(unsigned int i, SubVolume subvolume){
  this->subvolumes[i] = subvolume;
}

double Block::get_xmin() const {
  return (this->get_pos().x - this->get_dim().x/2.0);
}

double Block::get_xmax() const {
  return (this->get_pos().x + this->get_dim().x/2.0);
}

double Block::get_ymin() const {
  return (this->get_pos().y - this->get_dim().y/2.0);
}

double Block::get_ymax() const{
  return (this->get_pos().y + this->get_dim().y/2.0);
}

double Block::get_zmin() const{
  return (this->get_pos().z - this->get_dim().z/2.0);
}

double Block::get_zmax() const{
  return (this->get_pos().z + this->get_dim().z/2.0);
}

double Block::get_physical_length() const{
  return this->get_dim().z;
}
