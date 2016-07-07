#include <cmath>
#include <idmodel.h>

void DELTA::gen_delta(const Block& genblock, const unsigned int nr_periods, const double vertical_gap, const double horizontal_gap, const double block_separation){
  this->cs.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->ci.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->ce.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
  this->cd.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);

  Vector3D<double> dim = genblock.get_dim();

  this->cs.set_x(0.0);
  this->cs.set_z(+(vertical_gap + dim.z)/2.0);
  this->cs.set_ycenter(0.0);

  this->ci.set_x(0.0);
  this->ci.set_z(-(vertical_gap + dim.z)/2.0);
  this->ci.set_ycenter(0.0);

  this->cd.set_x(+(horizontal_gap + dim.x)/2.0);
  this->cd.set_z(0.0);
  this->cd.set_ycenter(0.0);

  this->ce.set_x(-(horizontal_gap + dim.x)/2.0);
  this->ce.set_z(0.0);
  this->ce.set_ycenter(0.0);

}

DELTA::DELTA(const Block& genblock, const unsigned int nr_periods, const double vertical_gap, const double horizontal_gap, const double block_separation){
  this->gen_delta(genblock, nr_periods, vertical_gap, horizontal_gap, block_separation);
};

Vector3D<double> DELTA::get_field(const Vector3D<double>& pos) const {
  Vector3D<double> field;
  field += this->cs.get_field(pos);
  field += this->ci.get_field(pos);
  field += this->cd.get_field(pos);
  field += this->ce.get_field(pos);

  Matrix3D<double> I = Matrix3D<double>::I();
  Matrix3D<double> rotz45 = I.set_rotation_z(-M_PI/4.0);
  field = rotz45*field;

  return field;
}

void DELTA::set_phase_cs(const double phase){
  this->cs.set_ycenter(phase);
}

void DELTA::set_phase_ci(const double phase){
  this->ci.set_ycenter(phase);
}
