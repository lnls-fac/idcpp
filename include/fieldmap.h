#ifndef FIELD_MAP_H
#define FIELD_MAP_H

#include <string>
#include <vector>
#include "vector3d.hpp"
#include "linalg.h"
#include "interpolation.h"

class FieldMap {

public:

	friend void  load_fieldmap(const std::string& fname, size_t& id, size_t& nx, double& x_min, double& x_max, size_t& nz, double& z_min, double& z_max);
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
	bool interpolation2d;
	alglib::spline3dinterpolant interpolant3d;
	alglib::spline2dinterpolant interpolant2d;

	FieldMap(const std::string& fname_, size_t id_ = 0);

	size_t           getid() const { return this->id; }
	Vector3D<double> field(const Vector3D<double>& pos) const;
	void             delete_data();

private:

	void read_fieldmap_from_file(const std::string& fname_);

};



#endif
