#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <set>
#include <API.h>
#include <linalg.h>
#include <interpolation.h>

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

	std::cout << "nr.points       : " << nr_points    << std::endl;
	std::cout << "nr.points.x_set : " << this->nx     << std::endl;
	std::cout << "nr.points.y_set : " << this->ny     << std::endl;
	std::cout << "nr.points.z_set : " << this->nz     << std::endl;
	std::cout << "min.x           : " << this->x_min  << std::endl;
	std::cout << "max.x           : " << this->x_max  << std::endl;
	std::cout << "min.y           : " << this->y_min  << std::endl;
	std::cout << "max.y           : " << this->y_max  << std::endl;
	std::cout << "min.z           : " << this->z_min  << std::endl;
	std::cout << "max.z           : " << this->z_max  << std::endl;
	std::cout << "dx              : " << this->dx     << std::endl;
	std::cout << "dy              : " << this->dy     << std::endl;
	std::cout << "dz              : " << this->dz     << std::endl;

	// throws exception in case dimensions do not agree
	if (nr_points != (this->nx * this->nz * this->ny)) {
		throw FieldMapException::inconsistent_dimensions;
	}

	this->map2d = (this->ny == 1) ? true : false;

	int data_array_size = 3*this->nx*this->ny*this->nz;
	alglib::real_1d_array x_array;
	alglib::real_1d_array y_array;
	alglib::real_1d_array z_array;
	alglib::real_1d_array data_array;
	x_array.setcontent(this->nx, &this->x_grid[0]);
	y_array.setcontent(this->ny, &this->y_grid[0]);
	z_array.setcontent(this->nz, &this->z_grid[0]);
	data_array.setcontent(data_array_size, &this->data[0]);

	if (this->map2d){
		std::cout << "2d" << std::endl;
		alglib::spline2dinterpolant interpolant2d;
		alglib::spline2dbuildbicubicv(x_array, nx, z_array, nz, data_array, 3, interpolant2d);
		this->interpolant2d = interpolant2d;
	} else{
		alglib::spline3dinterpolant interpolant3d;
		alglib::spline3dbuildtrilinearv(x_array, nx, y_array, ny, z_array, nz, data_array, 3, interpolant3d);
		this->interpolant3d = interpolant3d;
	}

	file.close();

}

Vector3D<double> FieldMap::field(const Vector3D<double>& pos) const{

	if (pos.x < this->x_min) { throw FieldMapException::out_of_range_x_min; }
	if (pos.x > this->x_max) { throw FieldMapException::out_of_range_x_max; }

	if (pos.z < this->z_min) { throw FieldMapException::out_of_range_z_min; }
	if (pos.z > this->z_max) { throw FieldMapException::out_of_range_z_max; }

	if (this->map2d){
		alglib::real_1d_array field2d;
		alglib::spline2dcalcv(interpolant2d, pos.x, pos.z, field2d);
		return Vector3D<double>(field2d[0], field2d[1], field2d[2]);

	} else {

		if (pos.y < this->y_min) { throw FieldMapException::out_of_range_y_min; }
		if (pos.y > this->y_max) { throw FieldMapException::out_of_range_y_max; }

		alglib::real_1d_array field3d;
		alglib::spline3dcalcv(interpolant3d, pos.x, pos.y, pos.z, field3d);
		return Vector3D<double>(field3d[0], field3d[1], field3d[2]);
	}
}
