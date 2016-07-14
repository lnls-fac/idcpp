#ifndef API_H
#define API_H

#include <string>
#include <vector>
#include <alglib/linalg.h>
#include <alglib/interpolation.h>
#include "vector3d.hpp"
#include "matrix3d.h"


struct InsertionDeviceException {
    enum type {
        success = 0,
        inconsistent_dimensions = 1,
        file_not_found = 2,
        invalid_magnetization = 3,
				invalid_shape = 4,
				field_symmetry_error = 5,
    };
};


class Subblock {
public:
  Subblock( Vector3D<double> dim_ = 0.0, Vector3D<double> pos_ = 0.0, double str_ = 1.0) : pos(pos_), dim(dim_), str(str_) {};

	Matrix3D<double> get_gmatrix( Vector3D<double> r);
  double           str;
  Vector3D<double> pos;
  Vector3D<double> dim;

};


class Block {
public:
	Block( Vector3D<double> mag_, Vector3D<double> dim_, Vector3D<double> pos_ = 0.0, double str_ = 1):
        mag(mag_) { subblocks.push_back(Subblock(dim_,pos_,str_)); }

	Block& add_subblock( Subblock sb) { subblocks.push_back(sb); }
  Block& add_subblock( Vector3D<double> dim_, Vector3D<double> pos_, double str_ = -1) { subblocks.push_back(Subblock(dim_,pos_,str_)); }
  Matrix3D<double>  get_gmatrix( Vector3D<double> r) ;
  Vector3D<double>  get_field( Vector3D<double> r) ;
  Vector3D<double>  get_mag() { return mag; }
  Vector3D<double>  get_pos() { return subblocks[0].pos; }
  Vector3D<double>  get_dim() { return subblocks[0].dim; }
  Subblock          get_subblock(unsigned int i) { return subblocks[i]; }
  void 						  set_mag(Vector3D<double> mag_);
  void  						set_pos(Vector3D<double> pos_);
  void  						set_dim(Vector3D<double> dim_);
  void              set_subblock(unsigned int i, Subblock subblock);

private:
  Vector3D<double>      mag;
  std::vector<Subblock> subblocks;

};


class Container {
public:
  Container() {}
  Container& add_element(Block& b) { blocks.push_back(b); return *this; }
  Vector3D<double> get_field( Vector3D<double>& r) ;
  std::vector<Vector3D<double> > get_field( std::vector<Vector3D<double> >& r) ;
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
	HalbachCassette(const HalbachCassette &obj);

	void gen_halbach_cassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing=0.0,  int N=4);
  HalbachCassette& 				set_x(double x);
  HalbachCassette& 				set_zcenter( double z=0.0);
  HalbachCassette& 				set_y(double y);
  Vector3D<double> 				get_pos();
  Vector3D<double> 				get_dim();
	Block&           				get_genblock() { return blocks[0];}
	unsigned int						get_number_of_periods() {return nr_periods;}
	double  								get_block_separation() {return spacing;}
	int	 	   								get_number_of_blocks_per_period() {return N;}

private:
	unsigned int nr_periods;
	double spacing;
	int N;
};

class EPU {
public:

  EPU() {};
  EPU(Block& genblock,  unsigned int nr_periods,  double magnetic_gap,  double cassette_separation=0.0,  double block_separation=0.0);
	EPU(const EPU &obj);

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
	unsigned int						get_number_of_periods()   {return nr_periods;}
	double									get_magnetic_gap()        {return magnetic_gap;}
	double  								get_cassette_separation() {return cassette_separation;}
	double  								get_block_separation()    {return block_separation;}

private:
	unsigned int nr_periods;
	double magnetic_gap;
	double cassette_separation;
	double block_separation;
};


class DELTA {
public:

  DELTA() {};
  DELTA(Block& genblock,  unsigned int nr_periods,  double vertical_gap,  double horizontal_gap,  double block_separation=0.0);
	DELTA(const DELTA &obj);

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
	unsigned int						get_number_of_periods()   {return nr_periods;}
	double									get_vertical_gap()        {return vertical_gap;}
	double  								get_horizontal_gap()      {return horizontal_gap;}
	double  								get_block_separation()    {return block_separation;}

private:
	unsigned int nr_periods;
	double vertical_gap;
	double horizontal_gap;
	double block_separation;
};


class Grid {
public:
	Grid() {};
	Grid(int nx, int ny, double x_max, double y_max);
  Grid(int nx, int ny, double x_min, double x_max, double y_min, double y_max);
	Grid(const Grid &obj);

	int nx;
  int ny;
	std::vector<double> x;
	std::vector<double> y;

private:
  double x_min, dx, x_max;
  double y_min, dy, y_max;

};


class Mask {
public:
	Mask() {};
  Mask(std::string shape, double width, double height);
  Mask(std::string filename);
	Mask(const Mask &obj);

	void load(std::string shape, double width, double height);
	void load(std::string filename);
  bool is_inside(Vector3D<> position) const;

private:

  double width;
  double height;
  std::vector<double> x;
  std::vector<double> y;
  std::string shape;
  alglib::spline1dinterpolant interpolant;

  void calc_interpolant();

};


class FieldMap {
public:

	FieldMap(const std::string& fname_, size_t id_ = 0);
	FieldMap(const FieldMap &obj);

	size_t id;
	size_t nx;
	size_t nz;
	double x_min, dx, x_max;
	double z_min, dz, z_max;
	double y;
	double physical_length;
	std::vector<double> x_grid;
	std::vector<double> z_grid;
	double *data;
	std::string fname;

	Vector3D<double> field(const Vector3D<double>& pos) const;
	std::vector<Vector3D<double> > field(const std::vector<Vector3D<double> >& pos) const;
	size_t  getid() const { return this->id; }
	void    delete_data();
	void    use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd);

private:

	alglib::spline2dinterpolant interpolant;
	void calc_interpolant();
	void read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions);

};


class FieldMap3D {
public:

	FieldMap3D() {};
	FieldMap3D(const std::string& fname_, size_t id_ = 0);
	FieldMap3D(const FieldMap3D &obj);

	size_t id;
	size_t nx;
	size_t ny;
	size_t nz;
	double x_min, dx, x_max;
	double y_min, dy, y_max;
	double z_min, dz, z_max;
	double physical_length;
	std::vector<double> x_grid;
	std::vector<double> y_grid;
	std::vector<double> z_grid;
	double *data;
	std::string fname;

	Vector3D<double> field(const Vector3D<double>& pos) const;
	std::vector<Vector3D<double> > field(const std::vector<Vector3D<double> >& pos) const;
	size_t           getid() const { return this->id; }
	void             delete_data();

private:

	std::vector<alglib::spline2dinterpolant> interpolants;
	void get_y_data(std::vector<double>& y_data, int y_index);
	void calc_interpolants();
	void read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions);
	Vector3D<double> field2D(const Vector3D<double>& pos, int interpolant_index) const;

};


class FieldMapContainer {
public:

	FieldMapContainer() {};
	FieldMapContainer(FieldMap fieldmap);
	FieldMapContainer(FieldMap3D fieldmap3D);
	FieldMapContainer(std::vector<FieldMap> fieldmaps);
	FieldMapContainer(std::string fieldmap_filename, bool fieldmap3D = false);
  FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool fieldmap3D = false);

	int nr_fieldmaps;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	double physical_length;

	Vector3D<double> field(const Vector3D<double>& pos) const;
	std::vector<Vector3D<double> > field(const std::vector<Vector3D<double> >& pos) const;
	void use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd);

private:

  std::vector<FieldMap> fieldmaps;
	FieldMap3D fieldmap3D;
	void set_attributes();

};


class KickMap {
public:

	KickMap() {};
	KickMap(double physical_length, std::vector<double> x, std::vector<double> y, std::vector<std::vector<double> > kick_x, std::vector<std::vector<double> > kicky);

	double physical_length;
	int nx;
	int ny;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

	void write_to_file(std::string filename);

};


class InsertionDevice{
public:

	typedef Vector3D<> (*FieldFunction)(Vector3D<>);

	InsertionDevice() {};
	InsertionDevice(FieldMap& fieldmap);
	InsertionDevice(std::vector<FieldMap>& fieldmaps);
	InsertionDevice(FieldMapContainer& fieldmap_container);
	InsertionDevice(EPU& epu);
	InsertionDevice(DELTA& delta);
	InsertionDevice(Vector3D<> (*function)(Vector3D<>));

	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	double physical_length;
	KickMap kickmap;

	Vector3D<double> field(Vector3D<double>& pos);
	std::vector<Vector3D<double> > field(std::vector<Vector3D<double> >& pos);
	void write_fieldmap_file(std::string filename, std::vector<double> x_vector, double y, std::vector<double> z_vector);
	void write_fieldmap_file(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
	void write_fieldmap_files(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
	void calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step);
	void calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step, std::string trajectories_filename);
	void write_kickmap_file(std::string filename);

private:
	FieldFunction field_function;
	FieldMapContainer fieldmaps;
	EPU epu;
	DELTA delta;
	int type;

};


void calc_brho(double energy, double& brho, double& beta);

void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds);

void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks, std::vector<std::vector<double> >& trajectory);
void runge_kutta(Vector3D<> (*function)(Vector3D<>), double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks, std::vector<std::vector<double> >& trajectory);

void _calc_kickmap(InsertionDevice& id, Grid grid, Mask mask, double energy, double rk_step, double z_min, double z_max, double physical_length, KickMap& kickmap, std::vector<std::vector<std::vector<double> > >& trajectories);
void _calc_kickmap(Vector3D<> (*function)(Vector3D<>), Grid grid, Mask mask, double energy, double rk_step, double z_min, double z_max, double physical_length, KickMap& kickmap, std::vector<std::vector<std::vector<double> > >& trajectories);

void _write_fieldmap_file(InsertionDevice& insertiondevice, std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
void _write_fieldmap_file(Vector3D<> (*function)(Vector3D<>), std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);

void _write_fieldmap_files(InsertionDevice& insertiondevice, std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
void _write_fieldmap_files(Vector3D<> (*function)(Vector3D<>), std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);

void _write_kickmap_file(std::string filename, double physical_length, std::vector<double>& x, std::vector<double>& y, std::vector<std::vector<double> >& kick_x, std::vector<std::vector<double> >& kick_y);

void save_trajectories(std::string filename, std::vector<std::vector<std::vector<double> > >& trajectories);

void save_trajectory(std::string filename, std::vector<std::vector<double> >& trajectory);

#endif
