#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <set>
#include <API.h>

FieldMap::FieldMap(const std::string& fname_, size_t id_) :
		id(id_),
		nx(0), nz(0), ny(0),
		x_min(0), dx(0), x_max(0),
		y_min(0), dy(0), y_max(0),
		z_min(0), dz(0), z_max(0),
		data(0)
{
	this->read_fieldmap_from_file(fname_);
}


void FieldMap::delete_data() {
	delete [] this->data;
}

inline
size_t FieldMap::ix(const double& x) const
{
	if (this->nx > 1){
		size_t ix = floor(0.5 + (x - this->x_min) / this->dx);
		if (ix < 0) { throw FieldMapException::out_of_range_x_min; }
		if (ix >= this->nx) { throw FieldMapException::out_of_range_x_max;  }
		return ix;
	} else return 0;
}

inline
size_t FieldMap::iy(const double& y) const
{
	if (this->ny > 1){
		size_t iy = floor(0.5 + (y - this->y_min) / this->dy);
		if (iy < 0) { throw FieldMapException::out_of_range_y_min;  }
		if (iy >= this->ny) { throw FieldMapException::out_of_range_y_max; }
		return iy;
	} else return 0;
}

inline
size_t FieldMap::iz(const double& z) const
{
	if (this->nz > 1){
		size_t iz = floor(0.5 + (z - this->z_min) / this->dz);
		if (iz < 0) {	throw FieldMapException::out_of_range_z_min; }
		if (iz >= this->nz) { throw FieldMapException::out_of_range_z_max;}
		return iz;
	} else return 0;
}

inline
double FieldMap::x(size_t ix) const
{
	if (ix < 0) { throw FieldMapException::out_of_range_x_min; }
	if (ix >= this->nx) { throw FieldMapException::out_of_range_x_max; }
	return (this->x_min + ix * this->dx);
}

inline
double FieldMap::y(size_t iy) const
{
	if (iy < 0) { throw FieldMapException::out_of_range_y_min; }
	if (iy >= this->ny) { throw FieldMapException::out_of_range_y_max; }
	return (this->y_min + iy * this->dy);
}

inline
double FieldMap::z(size_t iz) const
{
	if (iz < 0) { throw FieldMapException::out_of_range_z_min; }
	if (iz >= this->nz) { throw FieldMapException::out_of_range_z_max; }
	return (this->z_min + iz * this->dz);
}

void FieldMap::read_fieldmap_from_file(const std::string& fname_)
{

	const size_t capacity_inc = 10 * 1000;

	// opens file
	std::ifstream file(fname_.c_str());
	if (!file.is_open()) throw FieldMapException::file_not_found;

	// reads header section
	std::string word, line;
	file >> word; std::string nome_do_mapa; std::getline(file, nome_do_mapa);
	file >> word; std::string data_hora;    std::getline(file, data_hora);
	file >> word; std::string nome_arquivo; std::getline(file, nome_arquivo);
	file >> word; std::string nr_imas; 		std::getline(file, nr_imas);
	file >> word; std::string nome_ima;		std::getline(file, nome_ima);
	file >> word; std::string gap;  		std::getline(file, gap);
	file >> word; std::string gap_controle;	std::getline(file, gap_controle);
	file >> word; std::string comprimento;	std::getline(file, comprimento);
	file >> word; std::string corrente;		std::getline(file, corrente);
	file >> word; std::string centro_pos_z;	std::getline(file, centro_pos_z);
	file >> word; std::string centro_pos_x;	std::getline(file, centro_pos_x);
	file >> word; std::string rotacao;		std::getline(file, rotacao);
	std::getline(file, line);
	std::getline(file, line);
	std::getline(file, line);

	this->physical_length = std::atof(comprimento.c_str())/1000.0;

	// reads data section
	std::set<double> x_set, y_set, z_set;
	double x,y,z,bx,by,bz;
	size_t nr_points = 0;
	size_t capacity  = 0;
	if (this->data != 0) { delete [] this->data; } // deallocates memory if already set
	while (true) {
		if (nr_points >= capacity) {
			// reallocates more data space
			capacity += capacity_inc;
			this->data = (double*) std::realloc(this->data, capacity*3*sizeof(double));
		}
		file >> x >> y >> z >> bx >> by >> bz;
		data[3*nr_points+0] = bx;
		data[3*nr_points+1] = by;
		data[3*nr_points+2] = bz;
		if (file.eof()) break;
		x_set.insert(x/1000.0); y_set.insert(y/1000.0); z_set.insert(z/1000.0);
		nr_points++;
	}
	// updates object properties
	this->fname  = fname_;
	this->data   = (double*) std::realloc(this->data, nr_points*3*sizeof(double));
	this->nx     = x_set.size();
	this->ny     = y_set.size();
	this->nz     = z_set.size();
	this->x_min  = *(x_set.begin()) ;
	this->x_max  = *(x_set.rbegin()) ;
	this->y_min  = *(y_set.begin()) ;
	this->y_max  = *(y_set.rbegin()) ;
	this->z_min  = *(z_set.begin()) ;
	this->z_max  = *(z_set.rbegin()) ;
	this->dx     = (this->nx > 1) ? (this->x_max - this->x_min)/(this->nx - 1) : 0.0;
	this->dy     = (this->ny > 1) ? (this->y_max - this->y_min)/(this->ny - 1) : 0.0;
	this->dz     = (this->nz > 1) ? (this->z_max - this->z_min)/(this->nz - 1) : 0.0;
	this->x_grid.assign(x_set.begin(), x_set.end());
	this->y_grid.assign(y_set.begin(), y_set.end());
	this->z_grid.assign(z_set.begin(), z_set.end());

	// std::cout << "nr.points       : " << nr_points    << std::endl;
	// std::cout << "nr.points.x_set : " << this->nx     << std::endl;
	// std::cout << "nr.points.y_set : " << this->ny     << std::endl;
	// std::cout << "nr.points.z_set : " << this->nz     << std::endl;
	// std::cout << "min.x           : " << this->x_min  << std::endl;
	// std::cout << "max.x           : " << this->x_max  << std::endl;
	// std::cout << "min.y           : " << this->y_min  << std::endl;
	// std::cout << "max.y           : " << this->y_max  << std::endl;
	// std::cout << "min.z           : " << this->z_min  << std::endl;
	// std::cout << "max.z           : " << this->z_max  << std::endl;
	// std::cout << "dx              : " << this->dx     << std::endl;
	// std::cout << "dy              : " << this->dy     << std::endl;
	// std::cout << "dz              : " << this->dz     << std::endl;

	// throws exception in case dimensions do not agree
	if (nr_points != (this->nx * this->nz * this->ny)) {
		throw FieldMapException::inconsistent_dimensions;
	}

	file.close();

}

Vector3D<double> FieldMap::field2D(const Vector3D<double>& pos) const
{

	// calcs indices
	size_t ix1  = this->ix(pos.x);
	size_t iz1  = this->iz(pos.z);
	size_t ix2  = (ix1 == this->nx) ? ix1 : ix1 + 1;
	size_t iz2  = (iz1 == this->nz) ? iz1 : iz1 + 1;
	// gets field on grid points
	//      point1: (iz1, ix1), point2: (iz1, ix2), point3: (iz2, ix1), point4: (iz2, ix2)
	size_t p1 = iz1 * this->nx + ix1;
	size_t p2 = iz1 * this->nx + ix2;
	size_t p3 = iz2 * this->nx + ix1;
	size_t p4 = iz2 * this->nx + ix2;
	double x1 = this->x(ix1);
	double z1 = this->z(iz1);
	double x2  = this->x(ix2), z2 = this->z(iz2);
	double bx1 = this->data[3*p1+0], by1 = this->data[3*p1+1], bz1 = this->data[3*p1+2];
	double bx2 = this->data[3*p2+0], by2 = this->data[3*p2+1], bz2 = this->data[3*p2+2];
	double bx3 = this->data[3*p3+0], by3 = this->data[3*p3+1], bz3 = this->data[3*p3+2];
	double bx4 = this->data[3*p4+0], by4 = this->data[3*p4+1], bz4 = this->data[3*p4+2];
	// linear iterpolation: first along z (v_l=(p1+p3)/2, v_u=(p2+p4)/2), then along x (v = (v_l+v_u)/2)
	//
	//  (p2) ----- z ----- (p4)
	//   |                  |
	//   |                  |
	//   x                  x
	//   |                  |
	//   |                  |
	//  (p1) ----- z ----- (p3)
	//
	double bx,by,bz;
	if (iz1 == iz2) {
		if (ix1 == ix2) {
			bx = bx1;
			by = by1;
			bz = bz1;
		} else {
			bx = bx1 + (bx2 - bx1) * (pos.x - x1) / (x2 - x1);
			by = by1 + (by2 - by1) * (pos.x - x1) / (x2 - x1);
			bz = bz1 + (bz2 - bz1) * (pos.x - x1) / (x2 - x1);
		}
	} else {
		if (ix1 == ix2) {
			bx = bx1 + (bx3 - bx1) * (pos.z - z1) / (z2 - z1);
			by = by1 + (by3 - by1) * (pos.z - z1) / (z2 - z1);
			bz = bz1 + (bz3 - bz1) * (pos.z - z1) / (z2 - z1);
		} else {
			// BX
			double bx_l = bx1 + (bx3 - bx1) * (pos.z - z1) / (z2 - z1);
			double bx_u = bx2 + (bx4 - bx2) * (pos.z - z1) / (z2 - z1);
			bx   = bx_l + (bx_u - bx_l) * (pos.x - x1) / (x2 - x1);
			// BY
			double by_l = by1 + (by3 - by1) * (pos.z - z1) / (z2 - z1);
			double by_u = by2 + (by4 - by2) * (pos.z - z1) / (z2 - z1);
			by   = by_l + (by_u - by_l) * (pos.x - x1) / (x2 - x1);
			// BZ
			double bz_l = bz1 + (bz3 - bz1) * (pos.z - z1) / (z2 - z1);
			double bz_u = bz2 + (bz4 - bz2) * (pos.z - z1) / (z2 - z1);
			bz   = bz_l + (bz_u - bz_l) * (pos.x - x1) / (x2 - x1);
		}
	}

	return Vector3D<double>(bx,by,bz);

}

Vector3D<double> FieldMap::field3D(const Vector3D<double>& pos) const
{

	// calcs indices
	size_t ix1  = this->ix(pos.x);
	size_t iy1  = this->iy(pos.y);
	size_t iz1  = this->iz(pos.z);
	size_t ix2  = (ix1 == (this->nx - 1)) ? ix1 : ix1 + 1;
	size_t iy2  = (iy1 == (this->ny - 1)) ? iy1 : iy1 + 1;
	size_t iz2  = (iz1 == (this->nz - 1)) ? iz1 : iz1 + 1;

	// gets field on grid points
	//      point1: (ix1, iy1, iz1), point2: (ix2, iy1, iz1),
	//      point3: (ix1, iy2, iz1), point4: (ix2, iy2, iz1),
	//      point5: (ix1, iy1, iz2), point6: (ix2, iy1, iz2),
	//      point7: (ix1, iy2, iz2), point8: (ix2, iy2, iz2),
	size_t p1 = ix1 + iy1 * this->nx + iz1 * this->nx *this->ny;
	size_t p2 = ix2 + iy1 * this->nx + iz1 * this->nx *this->ny;
	size_t p3 = ix1 + iy2 * this->nx + iz1 * this->nx *this->ny;
	size_t p4 = ix2 + iy2 * this->nx + iz1 * this->nx *this->ny;
	size_t p5 = ix1 + iy1 * this->nx + iz2 * this->nx *this->ny;
	size_t p6 = ix2 + iy1 * this->nx + iz2 * this->nx *this->ny;
	size_t p7 = ix1 + iy2 * this->nx + iz2 * this->nx *this->ny;
	size_t p8 = ix2 + iy2 * this->nx + iz2 * this->nx *this->ny;
	double x1 = this->x(ix1);
	double y1 = this->y(iy1);
	double z1 = this->z(iz1);
	double x2 = this->x(ix2);
	double y2 = this->y(iy2);
	double z2 = this->z(iz2);

	double bx1 = this->data[3*p1+0], by1 = this->data[3*p1+1], bz1 = this->data[3*p1+2];
	double bx2 = this->data[3*p2+0], by2 = this->data[3*p2+1], bz2 = this->data[3*p2+2];
	double bx3 = this->data[3*p3+0], by3 = this->data[3*p3+1], bz3 = this->data[3*p3+2];
	double bx4 = this->data[3*p4+0], by4 = this->data[3*p4+1], bz4 = this->data[3*p4+2];
	double bx5 = this->data[3*p5+0], by5 = this->data[3*p5+1], bz5 = this->data[3*p5+2];
	double bx6 = this->data[3*p6+0], by6 = this->data[3*p6+1], bz6 = this->data[3*p6+2];
	double bx7 = this->data[3*p7+0], by7 = this->data[3*p7+1], bz7 = this->data[3*p7+2];
	double bx8 = this->data[3*p8+0], by8 = this->data[3*p8+1], bz8 = this->data[3*p8+2];

	double bx,by,bz;
	if (iz1 == iz2) {
		if (iy1 == iy2) {
			if (ix1 == ix2) {
				bx = bx1;
				by = by1;
				bz = bz1;

			} else {
				bx = bx1 + (bx2 - bx1) * (pos.x - x1) / (x2 - x1);
				by = by1 + (by2 - by1) * (pos.x - x1) / (x2 - x1);
				bz = bz1 + (bz2 - bz1) * (pos.x - x1) / (x2 - x1);

			}
		} else{
			if (ix1 == ix2){
				bx = bx1 + (bx3 - bx1) * (pos.y - y1) / (y2 - y1);
				by = by1 + (by3 - by1) * (pos.y - y1) / (y2 - y1);
				bz = bz1 + (bz3 - bz1) * (pos.y - y1) / (y2 - y1);

			} else{
				// BX
				double bx_l = bx1 + (bx3 - bx1) * (pos.y - y1) / (y2 - y1);
				double bx_u = bx2 + (bx4 - bx2) * (pos.y - y1) / (y2 - y1);
				bx   = bx_l + (bx_u - bx_l) * (pos.x - x1) / (x2 - x1);
				// BY
				double by_l = by1 + (by3 - by1) * (pos.y - y1) / (y2 - y1);
				double by_u = by2 + (by4 - by2) * (pos.y - y1) / (y2 - y1);
				by   = by_l + (by_u - by_l) * (pos.x - x1) / (x2 - x1);
				// BZ
				double bz_l = bz1 + (bz3 - bz1) * (pos.y - y1) / (y2 - y1);
				double bz_u = bz2 + (bz4 - bz2) * (pos.y - y1) / (y2 - y1);
				bz   = bz_l + (bz_u - bz_l) * (pos.x - x1) / (x2 - x1);

			}
		}
	} else{
		if (iy1 == iy2){
			if (ix1 == ix2){
				bx = bx1 + (bx5 - bx1) * (pos.z - z1) / (z2 - z1);
				by = by1 + (by5 - by1) * (pos.z - z1) / (z2 - z1);
				bz = bz1 + (bz5 - bz1) * (pos.z - z1) / (z2 - z1);

			} else{
				// BX
				double bx_l = bx1 + (bx5 - bx1) * (pos.z - z1) / (z2 - z1);
				double bx_u = bx2 + (bx6 - bx2) * (pos.z - z1) / (z2 - z1);
				bx   = bx_l + (bx_u - bx_l) * (pos.x - x1) / (x2 - x1);
				// BY
				double by_l = by1 + (by5 - by1) * (pos.z - z1) / (z2 - z1);
				double by_u = by2 + (by6 - by2) * (pos.z - z1) / (z2 - z1);
				by   = by_l + (by_u - by_l) * (pos.x - x1) / (x2 - x1);
				// BZ
				double bz_l = bz1 + (bz5 - bz1) * (pos.z - z1) / (z2 - z1);
				double bz_u = bz2 + (bz6 - bz2) * (pos.z - z1) / (z2 - z1);
				bz   = bz_l + (bz_u - bz_l) * (pos.x - x1) / (x2 - x1);

			}
		} else{
			if (ix1 == ix2){
				// BX
				double bx_l = bx1 + (bx5 - bx1) * (pos.z - z1) / (z2 - z1);
				double bx_u = bx3 + (bx7 - bx3) * (pos.z - z1) / (z2 - z1);
				bx   = bx_l + (bx_u - bx_l) * (pos.y - y1) / (y2 - y1);
				// BY
				double by_l = by1 + (by5 - by1) * (pos.z - z1) / (z2 - z1);
				double by_u = by3 + (by7 - by3) * (pos.z - z1) / (z2 - z1);
				by   = by_l + (by_u - by_l) * (pos.y - y1) / (y2 - y1);
				// BZ
				double bz_l = bz1 + (bz5 - bz1) * (pos.z - z1) / (z2 - z1);
				double bz_u = bz3 + (bz7 - bz3) * (pos.z - z1) / (z2 - z1);
				bz   = bz_l + (bz_u - bz_l) * (pos.y - y1) / (y2 - y1);

			} else{
				// y = y1
				// BX
				double bx_l_y1 = bx1 + (bx5 - bx1) * (pos.z - z1) / (z2 - z1);
				double bx_u_y1 = bx2 + (bx6 - bx2) * (pos.z - z1) / (z2 - z1);
				double bx_y1   = bx_l_y1 + (bx_u_y1 - bx_l_y1) * (pos.x - x1) / (x2 - x1);
				// BY
				double by_l_y1 = by1 + (by5 - by1) * (pos.z - z1) / (z2 - z1);
				double by_u_y1 = by2 + (by6 - by2) * (pos.z - z1) / (z2 - z1);
				double by_y1   = by_l_y1 + (by_u_y1 - by_l_y1) * (pos.x - x1) / (x2 - x1);
				// BZ
				double bz_l_y1 = bz1 + (bz5 - bz1) * (pos.z - z1) / (z2 - z1);
				double bz_u_y1 = bz2 + (bz6 - bz2) * (pos.z - z1) / (z2 - z1);
				double bz_y1   = bz_l_y1 + (bz_u_y1 - bz_l_y1) * (pos.x - x1) / (x2 - x1);

				// y = y2
				// BX
				double bx_l_y2 = bx3 + (bx7 - bx3) * (pos.z - z1) / (z2 - z1);
				double bx_u_y2 = bx4 + (bx8 - bx4) * (pos.z - z1) / (z2 - z1);
				double bx_y2   = bx_l_y2 + (bx_u_y2 - bx_l_y2) * (pos.x - x1) / (x2 - x1);
				// BY
				double by_l_y2 = by3 + (by7 - by3) * (pos.z - z1) / (z2 - z1);
				double by_u_y2 = by4 + (by8 - by4) * (pos.z - z1) / (z2 - z1);
				double by_y2   = by_l_y2 + (by_u_y2 - by_l_y2) * (pos.x - x1) / (x2 - x1);
				// BZ
				double bz_l_y2 = bz3 + (bz7 - bz3) * (pos.z - z1) / (z2 - z1);
				double bz_u_y2 = bz4 + (bz8 - bz4) * (pos.z - z1) / (z2 - z1);
				double bz_y2   = bz_l_y2 + (bz_u_y2 - bz_l_y2) * (pos.x - x1) / (x2 - x1);

				// Interpolation along y
				bx = bx_y1 + (bx_y2 - bx_y1) * (pos.y - y1) / (y2 - y1);
				by = by_y1 + (by_y2 - by_y1) * (pos.y - y1) / (y2 - y1);
				bz = bz_y1 + (bz_y2 - bz_y1) * (pos.y - y1) / (y2 - y1);

			}
		}
	}

	return Vector3D<double>(bx,by,bz);

}
