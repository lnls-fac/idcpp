#ifndef API_H
#define API_H

#include <string>
#include <vector>
#include "vector3d.hpp"
#include "linalg.h"
#include "interpolation.h"

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

class KickMap {

public:

	KickMap(FieldMap& fieldmap, Grid grid, Mask mask, double energy, double runge_kutta_step);
	KickMap(std::vector<FieldMap>& fieldmaps, Grid grid, Mask mask, double energy, double runge_kutta_step);
	KickMap(FieldMapContainer& fieldmap_container, Grid grid, Mask mask, double energy, double runge_kutta_step);

	void write_kickmap(std::string filename);

private:

	FieldMapContainer fieldmaps;
  Grid grid;
  Mask mask;
  double energy;
  double runge_kutta_step;

	std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

	void calc_kicks();

};

void runge_kutta(FieldMapContainer& fieldmaps, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks);


#endif
