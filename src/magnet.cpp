#include <api.h>

std::vector<Vector3D<double> > Magnet::field_vector(const std::vector<Vector3D<double> >& pos) const{
  std::vector<Vector3D<> > field;
  for (int i=0; i < pos.size(); i+=1){
    field.push_back(this->field(pos[i]));
  }
  return field;
}
