#include <vector>
#include <algorithm>
#include <halbachcassette.h>

void EPU::gen_epu(Block& genblock,  unsigned int nr_periods,  double magnetic_gap,  double cassette_separation,  double block_separation){
  this->csd.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->cse.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->cid.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
  this->cie.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);

  Vector3D<double> dim = genblock.get_dim();
  this->csd.set_x(+(cassette_separation + dim.x)/2.0);
  this->csd.set_z(+(magnetic_gap + dim.z)/2.0);
  this->csd.set_ycenter(0.0);

  this->cse.set_x(-(cassette_separation + dim.x)/2.0);
  this->cse.set_z(+(magnetic_gap + dim.z)/2.0);
  this->cse.set_ycenter(0.0);

  this->cie.set_x(-(cassette_separation + dim.x)/2.0);
  this->cie.set_z(-(magnetic_gap + dim.z)/2.0);
  this->cie.set_ycenter(0.0);

  this->cid.set_x(+(cassette_separation + dim.x)/2.0);
  this->cid.set_z(-(magnetic_gap + dim.z)/2.0);
  this->cid.set_ycenter(0.0);
}

EPU::EPU(Block& genblock,  unsigned int nr_periods,  double magnetic_gap,  double cassette_separation,  double block_separation){
  this->gen_epu(genblock, nr_periods, magnetic_gap, cassette_separation, block_separation);
};

Vector3D<double> EPU::field( Vector3D<double>& pos)  {
  Vector3D<double> field;
  field += this->csd.get_field(pos);
  field += this->cse.get_field(pos);
  field += this->cie.get_field(pos);
  field += this->cid.get_field(pos);
  return field;
}

void EPU::set_phase_csd( double phase){
  this->csd.set_ycenter(phase);
}

void EPU::set_phase_cie( double phase){
  this->cie.set_ycenter(phase);
}

double EPU::get_x_min(){
  std::vector<double> x_min_vector;
  x_min_vector.push_back(this->csd.get_pos().x - this->csd.get_dim().x/2.0);
  x_min_vector.push_back(this->cse.get_pos().x - this->cse.get_dim().x/2.0);
  x_min_vector.push_back(this->cie.get_pos().x - this->cie.get_dim().x/2.0);
  x_min_vector.push_back(this->cid.get_pos().x - this->cid.get_dim().x/2.0);
  double x_min = *std::min_element(x_min_vector.begin(), x_min_vector.end());
  return x_min;
}

double EPU::get_y_min(){
  std::vector<double> y_min_vector;
  y_min_vector.push_back(this->csd.get_pos().y - this->csd.get_dim().y/2.0);
  y_min_vector.push_back(this->cse.get_pos().y - this->cse.get_dim().y/2.0);
  y_min_vector.push_back(this->cie.get_pos().y - this->cie.get_dim().y/2.0);
  y_min_vector.push_back(this->cid.get_pos().y - this->cid.get_dim().y/2.0);
  double y_min = *std::min_element(y_min_vector.begin(), y_min_vector.end());
  return y_min;
}

double EPU::get_z_min(){
  std::vector<double> z_min_vector;
  z_min_vector.push_back(this->csd.get_pos().z - this->csd.get_dim().z/2.0);
  z_min_vector.push_back(this->cse.get_pos().z - this->cse.get_dim().z/2.0);
  z_min_vector.push_back(this->cie.get_pos().z - this->cie.get_dim().z/2.0);
  z_min_vector.push_back(this->cid.get_pos().z - this->cid.get_dim().z/2.0);
  double z_min = *std::min_element(z_min_vector.begin(), z_min_vector.end());
  return z_min;
}

double EPU::get_x_max(){
  std::vector<double> x_max_vector;
  x_max_vector.push_back(this->csd.get_pos().x + this->csd.get_dim().x/2.0);
  x_max_vector.push_back(this->cse.get_pos().x + this->cse.get_dim().x/2.0);
  x_max_vector.push_back(this->cie.get_pos().x + this->cie.get_dim().x/2.0);
  x_max_vector.push_back(this->cid.get_pos().x + this->cid.get_dim().x/2.0);
  double x_max = *std::min_element(x_max_vector.begin(), x_max_vector.end());
  return x_max;
}

double EPU::get_y_max(){
  std::vector<double> y_max_vector;
  y_max_vector.push_back(this->csd.get_pos().y + this->csd.get_dim().y/2.0);
  y_max_vector.push_back(this->cse.get_pos().y + this->cse.get_dim().y/2.0);
  y_max_vector.push_back(this->cie.get_pos().y + this->cie.get_dim().y/2.0);
  y_max_vector.push_back(this->cid.get_pos().y + this->cid.get_dim().y/2.0);
  double y_max = *std::min_element(y_max_vector.begin(), y_max_vector.end());
  return y_max;
}

double EPU::get_z_max(){
  std::vector<double> z_max_vector;
  z_max_vector.push_back(this->csd.get_pos().z + this->csd.get_dim().z/2.0);
  z_max_vector.push_back(this->cse.get_pos().z + this->cse.get_dim().z/2.0);
  z_max_vector.push_back(this->cie.get_pos().z + this->cie.get_dim().z/2.0);
  z_max_vector.push_back(this->cid.get_pos().z + this->cid.get_dim().z/2.0);
  double z_max = *std::min_element(z_max_vector.begin(), z_max_vector.end());
  return z_max;
}

double EPU::get_physical_length(){
  double physical_length = std::fabs(this->get_y_max() - this->get_y_min());
  return physical_length;
}
