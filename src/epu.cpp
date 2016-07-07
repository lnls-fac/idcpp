#include <idmodel.h>

void EPU::gen_epu(const Block& genblock, const unsigned int nr_periods, const double magnetic_gap, const double& cassette_separation, const double& block_separation){
  this->csd.gen_halbach_cassette(genblock,  Matrix3D<>::rotx90n, nr_periods, block_separation);
  this->cse.gen_halbach_cassette(genblock,  Matrix3D<>::rotx90n, nr_periods, block_separation);
  this->cid.gen_halbach_cassette(genblock,  Matrix3D<>::rotx90p, nr_periods, block_separation);
  this->cie.gen_halbach_cassette(genblock,  Matrix3D<>::rotx90p, nr_periods, block_separation);

  Vector3D<> dim = genblock.get_dim();
  this->csd.set_x(+(cassette_separation + dim[0])/2.0);
  this->csd.set_z(+(magnetic_gap + dim[2])/2.0);
  this->csd.set_ycenter(0);

  this->cse.set_x(-(cassette_separation + dim[0])/2.0);
  this->cse.set_z(+(magnetic_gap + dim[2])/2.0);
  this->cse.set_ycenter(0);

  this->cie.set_x(-(cassette_separation + dim[0])/2.0);
  this->cie.set_z(-(magnetic_gap + dim[2])/2.0);
  this->cie.set_ycenter(0);

  this->cid.set_x(+(cassette_separation + dim[0])/2.0);
  this->cid.set_z(-(magnetic_gap + dim[2])/2.0);
  this->cid.set_ycenter(0);
}

EPU::EPU(const Block& genblock, const unsigned int nr_periods, const double magnetic_gap, const double& cassette_separation, const double& block_separation){
  this->gen_epu(genblock, nr_periods, magnetic_gap, cassette_separation, block_separation);
};

template <typename T>
Vector3D<T> EPU::get_field(const Vector3D<T>& pos) const {
  Vector3D<T> field;
  field += this->csd.get_field(pos);
  field += this->cse.get_field(pos);
  field += this->cie.get_field(pos);
  field += this->cid.get_field(pos);
  return field;
}

void EPU::set_phase_csd(const double phase){
  this->csd.set_ycenter(phase);
}

void EPU::set_phase_cie(const double phase){
  this->cie.set_ycenter(phase);
}
