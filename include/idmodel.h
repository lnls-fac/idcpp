#ifndef ID_MODEL_H
#define ID_MODEL_H

#include <vector>
#include "vector3d.h"
#include "matrix3d.h"

class Subblock {

public:
  Subblock(const Vector3D<double>& dim_ = 0,
           const Vector3D<double>& pos_ = 0,
           const double& str_ = 1.0) :
           pos(pos_), dim(dim_), str(str_) {};
  Matrix3D<double> get_gmatrix(const Vector3D<double>& r) const;
  double          str;
  Vector3D<double> pos;
  Vector3D<double> dim;

};


class Block {
public:
  Block(const Vector3D<double>& mag_,
        const Vector3D<double>& dim_,
        const Vector3D<double>& pos_ = 0,
        const double& str_ = 1) :
        mag(mag_) { subblocks.push_back(Subblock(dim_,pos_,str_)); }
  Block& add_subblock(const Subblock& sb) { subblocks.push_back(sb); }
  Block& add_subblock(const Vector3D<double>& dim_,
                      const Vector3D<double>& pos_,
                      const double& str_ = -1) { subblocks.push_back(Subblock(dim_,pos_,str_)); }
  Matrix3D<double> get_gmatrix(const Vector3D<double>& t) const;
  Vector3D<double> get_field(const Vector3D<double>& t) const;
  Vector3D<double>&       set_mag()       { return mag; }
  const Vector3D<double>& get_mag() const { return mag; }
  Vector3D<double>&       set_pos()       { return subblocks[0].pos; }
  const Vector3D<double>& get_pos() const { return subblocks[0].pos; }
  Vector3D<double>&       set_dim()       { return subblocks[0].dim; }
  const Vector3D<double>& get_dim() const { return subblocks[0].dim; }
  const Subblock         set_subblock(unsigned int i)       { return subblocks[i]; }
  const Subblock         get_subblock(unsigned int i) const { return subblocks[i]; }
private:
  Vector3D<double>       mag;
  std::vector<Subblock> subblocks;
};

class Container {

public:
  Container() {}
  Container& add_element(const Block& b) { blocks.push_back(b); return *this; }
  Vector3D<double> get_field(const Vector3D<double>& r) const;
  std::vector<Vector3D<double> > get_field(const std::vector<Vector3D<double> >& rvec) const;
  int size() const { return blocks.size(); }
  const Block& get_item(int i) const { return blocks[i]; }
  Container& shift_pos(const Vector3D<double>& dr);
protected:
  std::vector<Block> blocks;

};

class HalbachCassette : public Container {
public:
  HalbachCassette() {};
  HalbachCassette(const Block& genblock, const Matrix3D<double>& rot, const unsigned int nr_periods, const double& spacing=0, const int& N=4);
  void gen_halbach_cassette(const Block& genblock, const Matrix3D<double>& rot, const unsigned int nr_periods, const double& spacing=0, const int& N=4);
  HalbachCassette& set_x(const double& x);
  HalbachCassette& set_ycenter(const double& y=0);
  HalbachCassette& set_z(const double& z);
};

class EPU {
public:
  EPU() {};
  EPU(const Block& genblock, const unsigned int nr_periods, const double magnetic_gap, const double& cassette_separation=0, const double& block_separation=0);
  void gen_epu(const Block& genblock, const unsigned int nr_periods, const double magnetic_gap, const double& cassette_separation, const double& block_separation);
  void set_phase_csd(const double phase);
  void set_phase_cie(const double phase);
  Vector3D<double> get_field(const Vector3D<double>& pos) const;
  HalbachCassette csd;
  HalbachCassette cse;
  HalbachCassette cid;
  HalbachCassette cie;
};

class DELTA {
public:
  DELTA() {};
  DELTA(const Block& genblock, const unsigned int nr_periods, const double vertical_gap, const double horizontal_gap, const double block_separation=0);
  void gen_delta(const Block& genblock, const unsigned int nr_periods, const double vertical_gap, const double horizontal_gap, const double block_separation=0);
  void set_phase_cs(const double phase);
  void set_phase_ci(const double phase);
  Vector3D<double> get_field(const Vector3D<double>& pos) const;
  HalbachCassette cs;
  HalbachCassette ci;
  HalbachCassette ce;
  HalbachCassette cd;
};


#endif
