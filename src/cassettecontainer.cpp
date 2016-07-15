#include <api.h>

CassetteContainer::CassetteContainer(HalbachCassette& cassette){
  this->cassettes.push_back(cassette);
}

CassetteContainer::CassetteContainer(std::vector<HalbachCassette> cassettes){
  this->cassettes = cassettes;
}

Vector3D<double> CassetteContainer::field( Vector3D<double>& r){
  Vector3D<double> field;
  for (int i=0; i < this->cassettes.size(); i+=1){
    field += this->cassettes[i].field(r);
  }
  return field;
}

std::vector<Vector3D<double> > CassetteContainer::field( std::vector<Vector3D<double> >& r){
  std::vector<Vector3D<> > field;
  for (int i=0; i < r.size(); i+=1){
    field.push_back(this->field(r[i]));
  }
  return field;
}

double CassetteContainer::get_x_min(){
  std::vector<double> x_min_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    x_min_vector.push_back(this->cassettes[i].get_first_block_pos().x - this->cassettes[i].get_block_dim().x/2.0);
  }
  double x_min = *std::min_element(x_min_vector.begin(), x_min_vector.end());
  return x_min;
}

double CassetteContainer::get_y_min(){
  std::vector<double> y_min_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    y_min_vector.push_back(this->cassettes[i].get_first_block_pos().y - this->cassettes[i].get_block_dim().y/2.0);
  }
  double y_min = *std::min_element(y_min_vector.begin(), y_min_vector.end());
  return y_min;
}

double CassetteContainer::get_z_min(){
  std::vector<double> z_min_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    z_min_vector.push_back(this->cassettes[i].get_first_block_pos().z - this->cassettes[i].get_block_dim().z/2.0);
  }
  double z_min = *std::min_element(z_min_vector.begin(), z_min_vector.end());
  return z_min;
}

double CassetteContainer::get_x_max(){
  std::vector<double> x_max_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    x_max_vector.push_back(this->cassettes[i].get_last_block_pos().x + this->cassettes[i].get_block_dim().x/2.0);
  }
  double x_max = *std::min_element(x_max_vector.begin(), x_max_vector.end());
  return x_max;
}

double CassetteContainer::get_y_max(){
  std::vector<double> y_max_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    y_max_vector.push_back(this->cassettes[i].get_last_block_pos().y + this->cassettes[i].get_block_dim().y/2.0);
  }
  double y_max = *std::min_element(y_max_vector.begin(), y_max_vector.end());
  return y_max;
}

double CassetteContainer::get_z_max(){
  std::vector<double> z_max_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    z_max_vector.push_back(this->cassettes[i].get_last_block_pos().z + this->cassettes[i].get_block_dim().z/2.0);
  }
  double z_max = *std::min_element(z_max_vector.begin(), z_max_vector.end());
  return z_max;
}

double CassetteContainer::get_physical_length(){
  double physical_length = std::fabs(this->get_z_max() - this->get_z_min());
  return physical_length;
}
