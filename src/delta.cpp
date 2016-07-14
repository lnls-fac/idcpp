#include <cmath>
#include <vector>
#include <algorithm>
#include <api.h>

void DELTA::gen_delta(Block& genblock,  unsigned int nr_periods,  double vertical_gap,  double horizontal_gap,  double block_separation){
  this->nr_periods = nr_periods;
  this->vertical_gap = vertical_gap;
  this->horizontal_gap = horizontal_gap;
  this->block_separation = block_separation;

  this->cs.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->cd.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  this->ci.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90p(), nr_periods, block_separation);
  this->ce.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90p(), nr_periods, block_separation);

  Vector3D<double> dim = genblock.get_dim();

  this->cs.set_x(0.0);
  this->cs.set_y(+(vertical_gap + dim.y)/2.0);
  this->cs.set_zcenter(0.0);

  this->ci.set_x(0.0);
  this->ci.set_y(-(vertical_gap + dim.y)/2.0);
  this->ci.set_zcenter(0.0);

  this->cd.set_x(+(horizontal_gap + dim.x)/2.0);
  this->cd.set_y(0.0);
  this->cd.set_zcenter(0.0);

  this->ce.set_x(-(horizontal_gap + dim.x)/2.0);
  this->ce.set_y(0.0);
  this->ce.set_zcenter(0.0);

}

DELTA::DELTA(Block& genblock,  unsigned int nr_periods,  double vertical_gap,  double horizontal_gap,  double block_separation){
  this->gen_delta(genblock, nr_periods, vertical_gap, horizontal_gap, block_separation);
};

DELTA::DELTA(const DELTA &obj){
  this->nr_periods = obj.nr_periods;
  this->vertical_gap = obj.vertical_gap;
  this->horizontal_gap = obj.horizontal_gap;
  this->block_separation = obj.block_separation;
  this->cs = obj.cs;
  this->ci = obj.ci;
  this->ce = obj.ce;
  this->cd = obj.cd;
}

Vector3D<double> DELTA::field( Vector3D<double>& pos)  {
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

void DELTA::set_phase_cs( double phase){
  this->cs.set_zcenter(phase);
}

void DELTA::set_phase_ci( double phase){
  this->ci.set_zcenter(phase);
}

double DELTA::get_x_min(){
  std::vector<double> x_min_vector;
  x_min_vector.push_back(this->cs.get_pos().x - this->cs.get_dim().x/2.0);
  x_min_vector.push_back(this->ci.get_pos().x - this->ci.get_dim().x/2.0);
  x_min_vector.push_back(this->cd.get_pos().x - this->cd.get_dim().x/2.0);
  x_min_vector.push_back(this->ce.get_pos().x - this->ce.get_dim().x/2.0);
  double x_min = *std::min_element(x_min_vector.begin(), x_min_vector.end());
  return x_min;
}

double DELTA::get_y_min(){
  std::vector<double> y_min_vector;
  y_min_vector.push_back(this->cs.get_pos().y - this->cs.get_dim().y/2.0);
  y_min_vector.push_back(this->ci.get_pos().y - this->ci.get_dim().y/2.0);
  y_min_vector.push_back(this->cd.get_pos().y - this->cd.get_dim().y/2.0);
  y_min_vector.push_back(this->ce.get_pos().y - this->ce.get_dim().y/2.0);
  double y_min = *std::min_element(y_min_vector.begin(), y_min_vector.end());
  return y_min;
}

double DELTA::get_z_min(){
  std::vector<double> z_min_vector;
  z_min_vector.push_back(this->cs.get_pos().z - this->cs.get_dim().z/2.0);
  z_min_vector.push_back(this->ci.get_pos().z - this->ci.get_dim().z/2.0);
  z_min_vector.push_back(this->cd.get_pos().z - this->cd.get_dim().z/2.0);
  z_min_vector.push_back(this->ce.get_pos().z - this->ce.get_dim().z/2.0);
  double z_min = *std::min_element(z_min_vector.begin(), z_min_vector.end());
  return z_min;
}

double DELTA::get_x_max(){
  std::vector<double> x_max_vector;
  x_max_vector.push_back(this->cs.get_pos().x + this->cs.get_dim().x/2.0);
  x_max_vector.push_back(this->ci.get_pos().x + this->ci.get_dim().x/2.0);
  x_max_vector.push_back(this->cd.get_pos().x + this->cd.get_dim().x/2.0);
  x_max_vector.push_back(this->ce.get_pos().x + this->ce.get_dim().x/2.0);
  double x_max = *std::min_element(x_max_vector.begin(), x_max_vector.end());
  return x_max;
}

double DELTA::get_y_max(){
  std::vector<double> y_max_vector;
  y_max_vector.push_back(this->cs.get_pos().y + this->cs.get_dim().y/2.0);
  y_max_vector.push_back(this->ci.get_pos().y + this->ci.get_dim().y/2.0);
  y_max_vector.push_back(this->cd.get_pos().y + this->cd.get_dim().y/2.0);
  y_max_vector.push_back(this->ce.get_pos().y + this->ce.get_dim().y/2.0);
  double y_max = *std::min_element(y_max_vector.begin(), y_max_vector.end());
  return y_max;
}

double DELTA::get_z_max(){
  std::vector<double> z_max_vector;
  z_max_vector.push_back(this->cs.get_pos().z + this->cs.get_dim().z/2.0);
  z_max_vector.push_back(this->ci.get_pos().z + this->ci.get_dim().z/2.0);
  z_max_vector.push_back(this->cd.get_pos().z + this->cd.get_dim().z/2.0);
  z_max_vector.push_back(this->ce.get_pos().z + this->ce.get_dim().z/2.0);
  double z_max = *std::min_element(z_max_vector.begin(), z_max_vector.end());
  return z_max;
}

double DELTA::get_physical_length(){
  double physical_length = std::fabs(this->get_z_max() - this->get_z_min());
  return physical_length;
}
