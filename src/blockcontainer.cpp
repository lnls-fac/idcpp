#include <vector>
#include <api.h>

Vector3D<double> BlockContainer::field( Vector3D<double>& r)  {
  Vector3D<double> f;
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    f += blocks[i].field(r);
  }
  return f;
}

std::vector<Vector3D<double> > BlockContainer::field( std::vector<Vector3D<double> >& r)  {
  std::vector<Vector3D<double> > f;
  for(unsigned int j = 0; j != r.size(); j++) {
    Vector3D<double> tf;
    for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
      tf += blocks[i].field(r[j]);
    }
    f.push_back(tf);
  }
  return f;
}

BlockContainer& BlockContainer::shift_pos( Vector3D<double> dr) {
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    Vector3D<double> pos = blocks[i].get_pos();
    blocks[i].set_pos(pos + dr);
  }
  return *this;
}
