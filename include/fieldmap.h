#ifndef FIELD_MAP_H
#define FIELD_MAP_H

#include <string>
#include <vector>
#include "vector3d.hpp"

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
	std::vector<double> x_grid;
	std::vector<double> y_grid;
	std::vector<double> z_grid;
	double *data;
	std::string fname;

	FieldMap(const std::string& fname_, size_t id_ = 0);

	size_t       getid() const { return this->id; }
	size_t       ix(const double& x) const;
	size_t       iy(const double& y) const;
	size_t       iz(const double& z) const;
	double       x(size_t ix) const;
	double       y(size_t iy) const;
	double       z(size_t iz) const;
	double 			 physical_length;
	Vector3D<double> pos(size_t ix, size_t iy) const;
	Vector3D<double> field2D(const Vector3D<double>& pos) const;
	Vector3D<double> field3D(const Vector3D<double>& pos) const;
	void         delete_data();

private:

	void read_fieldmap_from_file(const std::string& fname_);

};



#endif
