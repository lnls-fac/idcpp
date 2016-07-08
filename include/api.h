#ifndef API_H
#define API_H

#include <string>
#include <vector>
#include <alglib/linalg.h>
#include <alglib/interpolation.h>
#include "vector3d.hpp"
#include "idmodel.h"

struct GridException {
	enum type { inconsistent_dimensions = 0,};
};

struct MaskException {
  enum type { file_not_found = 0, undefined_shape = 1 };
};

struct FieldMapException {
	enum type { inconsistent_dimensions = 0, file_not_found = 1, };
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

	size_t           getid() const { return this->id; }
	void             delete_data();
	void 						 change_y_sign();

private:

	alglib::spline2dinterpolant interpolant;
	void calc_interpolant();
	void read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions);

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
	FieldMapContainer(std::vector<FieldMap> fieldmaps);
  FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool use_field_simmetry = false);

	Vector3D<double> field(const Vector3D<double>& pos) const;
	std::vector<Vector3D<double> > field(const std::vector<Vector3D<double> >& pos) const;

private:

  std::vector<FieldMap> fieldmaps;
	void set_attributes();

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

	InsertionDevice() {};
	InsertionDevice(FieldMap& fieldmap);
	InsertionDevice(std::vector<FieldMap>& fieldmaps);
	InsertionDevice(FieldMapContainer& fieldmap_container);
	InsertionDevice(EPU& epu);
	InsertionDevice(DELTA& delta);

	Vector3D<double> field(Vector3D<double>& pos);

	void write_fieldmap_file(std::string filename, std::vector<double> x_vector, double y, std::vector<double> z_vector);
	void write_fieldmap_files(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);

private:

	FieldMapContainer fieldmaps;
	EPU epu;
	DELTA delta;
	int type;

};

class KickMap {

public:

	KickMap(InsertionDevice& insertiondevice, Grid grid, Mask mask, double energy, double runge_kutta_step);
	void write_kickmap(std::string filename);

private:

	InsertionDevice insertiondevice;
  Grid grid;
  Mask mask;
  double energy;
  double runge_kutta_step;

	std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

	void calc_kicks();

};

void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks);

#endif
