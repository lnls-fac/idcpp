#ifndef API_H
#define API_H

#include <string>
#include <vector>
#include <alglib/linalg.h>
#include <alglib/interpolation.h>
#include "vector3d.hpp"
#include "halbachcassette.h"

struct GridException {
	enum type { inconsistent_dimensions = 0,};
};

struct MaskException {
  enum type { file_not_found = 0, undefined_shape = 1 };
};

struct FieldMapException {
	enum type { inconsistent_dimensions = 0, file_not_found = 1, field_symmetry = 2 };
};

class Grid {

public:

	int nx;
  int ny;
	std::vector<double> x;
	std::vector<double> y;

	Grid() {};
	Grid(int nx, int ny, double x_max, double y_max);
  Grid(int nx, int ny, double x_min, double x_max, double y_min, double y_max);
	Grid(const Grid &obj);

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

	FieldMap(const std::string& fname_, size_t id_ = 0);
	FieldMap(const FieldMap &obj);

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

	FieldMap3D() {};
	FieldMap3D(const std::string& fname_, size_t id_ = 0);
	FieldMap3D(const FieldMap3D &obj);

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

  int nr_fieldmaps;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	double physical_length;

	FieldMapContainer() {};
	FieldMapContainer(FieldMap fieldmap);
	FieldMapContainer(FieldMap3D fieldmap3D);
	FieldMapContainer(std::vector<FieldMap> fieldmaps);
	FieldMapContainer(std::string fieldmap_filename, bool fieldmap3D = false);
  FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool fieldmap3D = false);

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
	void write_to_file(std::string filename);

private:

	double physical_length;
	int nx;
	int ny;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

};

class InsertionDevice{

public:

	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	double physical_length;
	KickMap kickmap;
	
	InsertionDevice() {};
	InsertionDevice(FieldMap& fieldmap);
	InsertionDevice(std::vector<FieldMap>& fieldmaps);
	InsertionDevice(FieldMapContainer& fieldmap_container);
	InsertionDevice(EPU& epu);
	InsertionDevice(DELTA& delta);

	Vector3D<double> field(Vector3D<double>& pos);
	std::vector<Vector3D<double> > field(std::vector<Vector3D<double> >& pos);

	void write_fieldmap_file(std::string filename, std::vector<double> x_vector, double y, std::vector<double> z_vector);
	void write_fieldmap_files(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
	void write_fieldmap3D_file(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
	void calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step);
	void write_kickmap_file(std::string filename);

private:

	FieldMapContainer fieldmaps;
	EPU epu;
	DELTA delta;
	int type;

};

void calc_brho(double energy, double& brho, double& beta);
void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds);
void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks);

#endif
