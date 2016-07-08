#ifndef ID_MODEL_H
#define ID_MODEL_H

#include <vector>
#include "vector3d.hpp"
#include "matrix3d.h"

class Subblock {

public:
  Subblock( Vector3D<double> dim_ = 0.0,
            Vector3D<double> pos_ = 0.0,
            double str_ = 1.0) :
           pos(pos_), dim(dim_), str(str_) {};
  Matrix3D<double> get_gmatrix( Vector3D<double> r) ;
  double          str;
  Vector3D<double> pos;
  Vector3D<double> dim;

};


class Block {
public:
  Block( Vector3D<double> mag_,
         Vector3D<double> dim_,
         Vector3D<double> pos_ = 0.0,
         double str_ = 1) :
        mag(mag_) { subblocks.push_back(Subblock(dim_,pos_,str_)); }
  Block& add_subblock( Subblock sb) { subblocks.push_back(sb); }
  Block& add_subblock( Vector3D<double> dim_,
                       Vector3D<double> pos_,
                       double str_ = -1) { subblocks.push_back(Subblock(dim_,pos_,str_)); }
  Matrix3D<double> get_gmatrix( Vector3D<double> t) ;
  Vector3D<double> get_field( Vector3D<double> t) ;
  Vector3D<double>&      set_mag() { return mag; }
  Vector3D<double>&      get_mag() { return mag; }
  Vector3D<double>&      set_pos() { return subblocks[0].pos; }
  Vector3D<double>&      get_pos() { return subblocks[0].pos; }
  Vector3D<double>&      set_dim() { return subblocks[0].dim; }
  Vector3D<double>&      get_dim() { return subblocks[0].dim; }
  Subblock               set_subblock(unsigned int i)       { return subblocks[i]; }
  Subblock               get_subblock(unsigned int i)  { return subblocks[i]; }
private:
  Vector3D<double>       mag;
  std::vector<Subblock> subblocks;
};

class Container {

public:
  Container() {}
  Container& add_element(Block& b) { blocks.push_back(b); return *this; }
  Vector3D<double> get_field( Vector3D<double>& r) ;
  std::vector<Vector3D<double> > get_field( std::vector<Vector3D<double> >& rvec) ;
  int size() { return blocks.size(); }
  Block& get_item(int i) { return blocks[i]; }
  Container& shift_pos( Vector3D<double> dr);
protected:
  std::vector<Block> blocks;

};

class HalbachCassette : public Container {
public:
  HalbachCassette() {};
  HalbachCassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing=0.0,  int N=4);
  void gen_halbach_cassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing=0.0,  int N=4);
  HalbachCassette& set_x(double x);
  HalbachCassette& set_ycenter( double y=0.0);
  HalbachCassette& set_z(double z);
  Vector3D<double>& get_pos();
  Vector3D<double>& get_dim();
};

class EPU {
public:

  EPU() {};
  EPU(Block& genblock,  unsigned int nr_periods,  double magnetic_gap,  double cassette_separation=0.0,  double block_separation=0.0);
  void gen_epu(Block& genblock,  unsigned int nr_periods,  double magnetic_gap,  double cassette_separation,  double block_separation);
  void set_phase_csd( double phase);
  void set_phase_cie( double phase);
  Vector3D<double> field( Vector3D<double>& pos) ;
  HalbachCassette csd;
  HalbachCassette cse;
  HalbachCassette cid;
  HalbachCassette cie;
  double get_x_min();
  double get_x_max();
  double get_y_min();
  double get_y_max();
  double get_z_min();
  double get_z_max();
  double get_physical_length();
};

class DELTA {
public:

  DELTA() {};
  DELTA(Block& genblock,  unsigned int nr_periods,  double vertical_gap,  double horizontal_gap,  double block_separation=0.0);
  void gen_delta(Block& genblock,  unsigned int nr_periods,  double vertical_gap,  double horizontal_gap,  double block_separation);
  void set_phase_cs( double phase);
  void set_phase_ci( double phase);
  Vector3D<double> field( Vector3D<double>& pos) ;
  HalbachCassette cs;
  HalbachCassette ci;
  HalbachCassette ce;
  HalbachCassette cd;
  double get_x_min();
  double get_x_max();
  double get_y_min();
  double get_y_max();
  double get_z_min();
  double get_z_max();
  double get_physical_length();
};

#endif
