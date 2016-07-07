#include <vector>
#include <idmodel.h>

template <typename T>
Vector3D<T> Container::get_field(const Vector3D<T>& r) const {
  Vector3D<T> f;
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    f += blocks[i].get_field<T>(r);
  }
  return f;
}

template <typename T>
std::vector<Vector3D<T>> Container::get_field(const std::vector<Vector3D<T>>& rvec) const {
  std::vector<Vector3D<T>> f;
  for(unsigned int j = 0; j != rvec.size(); j++) {
    Vector3D<T> tf;
    for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
      tf += blocks[i].get_field<T>(rvec[j]);
    }
    f.push_back(tf);
  }
  return f;
}

Container& Container::shift_pos(const Vector3D<double>& dr) {
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos() = blocks[i].set_pos() + dr;
  }
  return *this;
}


template <unsigned int N>
std::ostream& operator <<(std::ostream& out, const HalbachCassette<N>& v) {

  char buffer[256];
  auto n = v.size();
  std::string str;
  str = ""; std::sprintf(buffer, "posx:%+.3e",v[0].get_pos()[0]); str += buffer;
  out << str << std::endl;
  str = ""; std::sprintf(buffer, "posz:%+.3e",v[0].get_pos()[2]); str += buffer;
  out << str << std::endl;
  out << " nr |   posy   |   magx       magy       magz   " << std::endl;
  out << "================================================" << std::endl;
  for(auto i=0; i<n-1; ++i) {
    str = "";
    std::sprintf(buffer, "%04i|", i+1); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e|", v[i].get_pos()[1]); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v[i].get_mag()[0]); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v[i].get_mag()[1]); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e",  v[i].get_mag()[2]); str += std::string(buffer);
    out << str << std::endl;
  }
  //std::string str;
  str = "";
  std::sprintf(buffer, "%04i|", n); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e|", v[n-1].get_pos()[1]); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v[n-1].get_mag()[0]); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v[n-1].get_mag()[1]); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e",  v[n-1].get_mag()[2]); str += std::string(buffer);
  out << str;
  return out;
}

template <unsigned int N>
void HalbachCassette<N>::gen_halbach_cassette(const Block& genblock, const Matrix3D<double>& rot, const unsigned int nr_periods, const double& spacing) {
  const double dy = spacing + genblock.get_dim()[1];
  Block block = genblock;
  this->blocks.clear();
  for(unsigned int i=0; i<nr_periods; ++i) {
    block.set_mag() = genblock.get_mag();
    for(unsigned int j=0; j<N; ++j) {
      this->blocks.push_back(block);           // adds block to container
      block.set_pos()[1] += dy;                // shifts longitudinally
      block.set_mag() = rot * block.get_mag(); // rotates magetization vector
    }
  }
};

template <unsigned int N>
HalbachCassette<N>::HalbachCassette(const Block& genblock, const Matrix3D<double>& rot, const unsigned int nr_periods, const double& spacing) {
  this->gen_halbach_cassette(genblock, rot, nr_periods, spacing);
};

template <unsigned int N>
HalbachCassette<N>& HalbachCassette<N>::set_x(const double& x) {
  for(auto i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos()[0] = x;
  }
  return *this;
}

template <unsigned int N>
HalbachCassette<N>& HalbachCassette<N>::set_ycenter(const double& y) {
  double ym = 0.5 * (this->blocks.front().get_pos()[1] + this->blocks.back().get_pos()[1]);
  this->shift_pos(Vector3D<>(0,y - ym,0));
  return *this;
}

template <unsigned int N>
HalbachCassette<N>& HalbachCassette<N>::set_z(const double& z) {
  for(auto i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos()[2] = z;
  }
  return *this;
}
