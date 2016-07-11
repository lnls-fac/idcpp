#include <vector>
#include <halbachcassette.h>

Vector3D<double> Container::get_field( Vector3D<double>& r)  {
  Vector3D<double> f;
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    f += blocks[i].get_field(r);
  }
  return f;
}

std::vector<Vector3D<double> > Container::get_field( std::vector<Vector3D<double> >& rvec)  {
  std::vector<Vector3D<double> > f;
  for(unsigned int j = 0; j != rvec.size(); j++) {
    Vector3D<double> tf;
    for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
      tf += blocks[i].get_field(rvec[j]);
    }
    f.push_back(tf);
  }
  return f;
}

Container& Container::shift_pos( Vector3D<double> dr) {
  for(std::vector<Block>::size_type i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos() += dr;
  }
  return *this;
}


std::ostream& operator <<(std::ostream& out, HalbachCassette& v) {

  char buffer[256];
  auto n = v.size();
  std::string str;
  str = ""; std::sprintf(buffer, "posx:%+.3e",v.get_item(0).get_pos().x); str += buffer;
  out << str << std::endl;
  str = ""; std::sprintf(buffer, "posz:%+.3e",v.get_item(0).get_pos().z); str += buffer;
  out << str << std::endl;
  out << " nr |   posy   |   magx       magy       magz   " << std::endl;
  out << "================================================" << std::endl;
  for(auto i=0; i<n-1; ++i) {
    str = "";
    std::sprintf(buffer, "%04i|", i+1); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e|", v.get_item(i).get_pos().y); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v.get_item(i).get_mag().x); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e ", v.get_item(i).get_mag().y); str += std::string(buffer);
    std::sprintf(buffer, "%+.3e",  v.get_item(i).get_mag().z); str += std::string(buffer);
    out << str << std::endl;
  }
  //std::string str;
  str = "";
  std::sprintf(buffer, "%04i|", n); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e|", v.get_item(n-1).get_pos().y); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v.get_item(n-1).get_mag().x); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e ", v.get_item(n-1).get_mag().y); str += std::string(buffer);
  std::sprintf(buffer, "%+.3e",  v.get_item(n-1).get_mag().z); str += std::string(buffer);
  out << str;
  return out;
}

void HalbachCassette::gen_halbach_cassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing,  int N) {
  double dy = spacing + genblock.get_dim().y;
  Block block = genblock;
  this->blocks.clear();
  for(unsigned int i=0; i<nr_periods; ++i) {
    block.set_mag() = genblock.get_mag();
    for(unsigned int j=0; j<N; ++j) {
      this->blocks.push_back(block);           // adds block to container
      block.set_pos().y += dy;                // shifts longitudinally
      block.set_mag() = rot * block.get_mag(); // rotates magetization vector
    }
  }
};

HalbachCassette::HalbachCassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing,  int N) {
  this->gen_halbach_cassette(genblock, rot, nr_periods, spacing, N);
};

HalbachCassette& HalbachCassette::set_x( double x) {
  for(auto i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos().x = x;
  }
  return *this;
}

HalbachCassette& HalbachCassette::set_ycenter( double y) {
  double ym = 0.5 * (this->blocks.front().get_pos().y + this->blocks.back().get_pos().y);
  this->shift_pos(Vector3D<double>(0,y - ym,0));
  return *this;
}

HalbachCassette& HalbachCassette::set_z( double z) {
  for(auto i = 0; i != blocks.size(); i++) {
    blocks[i].set_pos().z = z;
  }
  return *this;
}

Vector3D<double>& HalbachCassette::get_pos(){
  Vector3D<double>& pos = this->blocks[0].get_pos();
  double ym = 0.5 * (this->blocks.front().get_pos().y + this->blocks.back().get_pos().y);
  pos.y = ym;
  return pos;
}

Vector3D<double>& HalbachCassette::get_dim(){
  Vector3D<double>& dim = this->blocks[0].get_dim();
  double ylength = std::fabs(this->blocks.front().get_pos().y - this->blocks.back().get_pos().y);
  dim.y = ylength;
  return dim;
}
