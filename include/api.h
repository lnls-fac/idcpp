#ifndef API_H
#define API_H

#include <vector>
#include "vector3d.hpp"
#include "linalg.h"
#include "interpolation.h"
#include "fieldmap.h"

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

class KickMap {

public:

	KickMap(std::vector<FieldMap> fieldmaps, Grid grid, Mask mask, double energy, double runge_kutta_step);

	void write_kickmap(std::string filename);

private:

	std::vector<FieldMap> fieldmaps;
  Grid grid;
  Mask mask;
  double energy;
  double runge_kutta_step;

	std::vector<std::vector<double>> kick_x;
	std::vector<std::vector<double>> kick_y;

	void calc_kicks();

};

void runge_kutta(std::vector<FieldMap>& fieldmaps, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks);


#endif
