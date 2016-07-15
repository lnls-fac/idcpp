#include <cmath>
#include <vector>
#include <api.h>

Subblock::Subblock(const Subblock &obj){
  this->dim = obj.dim;
  this->pos = obj.pos;
  this->str = obj.str;
}

Matrix3D<double> Subblock::get_gmatrix( Vector3D<double> r)  {
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
  this->subblocks = obj.subblocks;
}

Matrix3D<double> Block::get_gmatrix( Vector3D<double> r)  {
  Matrix3D<double> m;
  for(std::vector<Subblock>::size_type i = 0; i != subblocks.size(); i++) {
    m = subblocks[i].get_gmatrix(r);
  }
  return m;
}


Vector3D<double> Block::field( Vector3D<double> r)  {
  return this->get_gmatrix(r) * this->mag;
}

void Block::set_mag(Vector3D<double> mag_){
  this->mag = mag_;
}

void Block::set_pos(Vector3D<double> pos_){
  this->subblocks[0].pos = pos_;
}

void Block::set_dim(Vector3D<double> dim_){
  this->subblocks[0].dim = dim_;
}

void Block::set_subblock(unsigned int i, Subblock subblock){
  this->subblocks[i] = subblock;
}
