#include <vector>
#include <api.h>

std::ostream& operator <<(std::ostream& out, HalbachCassette& v) {

  char buffer[256];
  auto n = v.size();
  std::string str;
  str = ""; std::sprintf(buffer, "posx:%+.3e",v.get_item(0).get_pos().x); str += buffer;
  out << str << std::endl;
  str = ""; std::sprintf(buffer, "posy:%+.3e",v.get_item(0).get_pos().y); str += buffer;
  out << str << std::endl;
  out << " nr |   posy   |   magx       magy       magz   " << std::endl;
  out << "================================================" << std::endl;
  for(auto i=0; i<n-1; ++i) {
    str = "";
    std::sprintf(buffer, "%04i|", i+1); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e|", v.get_item(i).get_pos().z); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v.get_item(i).get_mag().x); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v.get_item(i).get_mag().z); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e",  v.get_item(i).get_mag().y); str += std::string(buffer);
    out << str << std::endl;
  }
  //std::string str;
  str = "";
  std::sprintf(buffer, "%04i|", n); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e|", v.get_item(n-1).get_pos().z); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v.get_item(n-1).get_mag().x); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v.get_item(n-1).get_mag().z); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e",  v.get_item(n-1).get_mag().y); str += std::string(buffer);
  out << str;
  return out;
}

void HalbachCassette::gen_halbach_cassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing,  int N) {
  this->nr_periods = nr_periods; this->spacing = spacing; this->N = N;

  double dz = spacing + genblock.get_dim().z;
  Block block = genblock;
  this->blocks.clear();

  for(unsigned int i=0; i<nr_periods; ++i) {
    block.set_mag(genblock.get_mag());
    for(unsigned int j=0; j<N; ++j) {
      this->blocks.push_back(block);                                           // adds block to block container
      Vector3D<double> pos = block.get_pos(); pos.z += dz; block.set_pos(pos); // shifts longitudinally
      block.set_mag(rot * block.get_mag());                                    // rotates magetization vector
    }
  }

};

HalbachCassette::HalbachCassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing,  int N) {
  this->gen_halbach_cassette(genblock, rot, nr_periods, spacing, N);
};

HalbachCassette::HalbachCassette(const HalbachCassette &obj){
  this->blocks      = obj.blocks;
  this->nr_periods  = obj.nr_periods;
  this->spacing     = obj.spacing;
  this->N           = obj.N;
}

HalbachCassette& HalbachCassette::set_xcenter( double x) {
  for(auto i = 0; i != blocks.size(); i++) {
    Vector3D<double> pos = blocks[i].get_pos(); pos.x = x; blocks[i].set_pos(pos);
  }
  return *this;
}

HalbachCassette& HalbachCassette::set_zcenter( double z) {
  double zm = 0.5 * (this->blocks.front().get_pos().z + this->blocks.back().get_pos().z);
  this->shift_pos(Vector3D<double>(0.0, 0.0, z - zm));
  return *this;
}

HalbachCassette& HalbachCassette::set_ycenter( double y) {
  for(auto i = 0; i != blocks.size(); i++) {
    Vector3D<double> pos = blocks[i].get_pos(); pos.y = y; blocks[i].set_pos(pos);
  }
  return *this;
}

HalbachCassette& HalbachCassette::set_center_pos(Vector3D<double> pos){
  Vector3D<double> old_pos = this->get_center_pos();
  Vector3D<double> dr = pos - old_pos;
  this->shift_pos(dr);
  return *this;
}

Vector3D<double> HalbachCassette::get_center_pos() const{
  Vector3D<double> center = this->blocks[0].get_pos();
  double zm = 0.5 * (this->blocks.front().get_pos().z + this->blocks.back().get_pos().z);
  center.z = zm;
  return center;
}

Vector3D<double> HalbachCassette::get_dim() const{
  Vector3D<double> dim = this->blocks[0].get_dim();
  double zlength = std::fabs(this->blocks.front().get_pos().z - this->blocks.back().get_pos().z);
  dim.z = zlength;
  return dim;
}
