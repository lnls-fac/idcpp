#include <api.h>

CassetteContainer::CassetteContainer(HalbachCassette& cassette){
  this->cassettes.push_back(cassette);
}

CassetteContainer::CassetteContainer(std::vector<HalbachCassette> cassettes){
  this->cassettes = cassettes;
}

Vector3D<double> CassetteContainer::field(const Vector3D<double>& r) const {
  Vector3D<double> field;
  for (int i=0; i < this->cassettes.size(); i+=1){
    field += this->cassettes[i].field(r);
  }
  return field;
}

double CassetteContainer::get_xmin() const{
  std::vector<double> xmin_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    xmin_vector.push_back(this->cassettes[i].get_xmin());
  }
  double xmin = *std::min_element(xmin_vector.begin(), xmin_vector.end());
  return xmin;
}

double CassetteContainer::get_xmax() const{
  std::vector<double> xmax_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    xmax_vector.push_back(this->cassettes[i].get_xmax());
  }
  double xmax = *std::max_element(xmax_vector.begin(), xmax_vector.end());
  return xmax;
}

double CassetteContainer::get_ymin() const{
  std::vector<double> ymin_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    ymin_vector.push_back(this->cassettes[i].get_ymin());
  }
  double ymin = *std::min_element(ymin_vector.begin(), ymin_vector.end());
  return ymin;
}

double CassetteContainer::get_ymax() const{
  std::vector<double> ymax_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    ymax_vector.push_back(this->cassettes[i].get_ymax());
  }
  double ymax = *std::max_element(ymax_vector.begin(), ymax_vector.end());
  return ymax;
}

double CassetteContainer::get_zmin() const{
  std::vector<double> zmin_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    zmin_vector.push_back(this->cassettes[i].get_zmin());
  }
  double zmin = *std::min_element(zmin_vector.begin(), zmin_vector.end());
  return zmin;
}

double CassetteContainer::get_zmax() const{
  std::vector<double> zmax_vector;
  for (int i=0; i < this->cassettes.size(); i+=1){
    zmax_vector.push_back(this->cassettes[i].get_zmax());
  }
  double zmax = *std::max_element(zmax_vector.begin(), zmax_vector.end());
  return zmax;
}

double CassetteContainer::get_physical_length() const{
  double physical_length = std::fabs(this->get_zmax() - this->get_zmin());
  return physical_length;
}
