#include <api.h>

Vector3D<double> BlockContainer::field(const Vector3D<double>& r) const {
  Vector3D<double> f;
  for(std::vector<Block>::size_type i = 0; i != this->blocks.size(); i++) {
    f += this->blocks[i].field(r);
  }
  return f;
}

BlockContainer& BlockContainer::shift_pos( Vector3D<double> dr) {
  for(std::vector<Block>::size_type i = 0; i != this->blocks.size(); i++) {
    Vector3D<double> pos = this->blocks[i].get_pos();
    this->blocks[i].set_pos(pos + dr);
  }
  return *this;
}

double BlockContainer::get_xmin() const{
  std::vector<double> xmin_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    xmin_vector.push_back(this->blocks[i].get_xmin());
  }
  double xmin = *std::min_element(xmin_vector.begin(), xmin_vector.end());
  return xmin;
}

double BlockContainer::get_xmax() const{
  std::vector<double> xmax_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    xmax_vector.push_back(this->blocks[i].get_xmax());
  }
  double xmax = *std::max_element(xmax_vector.begin(), xmax_vector.end());
  return xmax;
}

double BlockContainer::get_ymin() const{
  std::vector<double> ymin_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    ymin_vector.push_back(this->blocks[i].get_ymin());
  }
  double ymin = *std::min_element(ymin_vector.begin(), ymin_vector.end());
  return ymin;
}

double BlockContainer::get_ymax() const{
  std::vector<double> ymax_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    ymax_vector.push_back(this->blocks[i].get_ymax());
  }
  double ymax = *std::max_element(ymax_vector.begin(), ymax_vector.end());
  return ymax;
}

double BlockContainer::get_zmin() const{
  std::vector<double> zmin_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    zmin_vector.push_back(this->blocks[i].get_zmin());
  }
  double zmin = *std::min_element(zmin_vector.begin(), zmin_vector.end());
  return zmin;
}

double BlockContainer::get_zmax() const{
  std::vector<double> zmax_vector;
  for (int i=0; i < this->blocks.size(); i+=1){
    zmax_vector.push_back(this->blocks[i].get_zmax());
  }
  double zmax = *std::max_element(zmax_vector.begin(), zmax_vector.end());
  return zmax;
}

double BlockContainer::get_physical_length() const{
  double physical_length = std::fabs(this->get_zmax() - this->get_zmin());
  return physical_length;
}
