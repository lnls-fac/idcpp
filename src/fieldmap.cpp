#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <set>
#include <linalg.h>
#include <interpolation.h>
#include <stdlib.h>
#include <vector3d.hpp>
#include <fieldmap.h>

struct FieldMapException {
	enum type {
		success                       = 0,
		inconsistent_dimensions       = 1,
		out_of_range_x_min            = 2,
		out_of_range_x_max            = 3,
		out_of_range_z_min            = 4,
		out_of_range_z_max            = 5,
		file_not_found                = 6,
		out_of_range_y_min            = 7,
		out_of_range_y_max            = 8
	};
};

typedef std::vector<double> state_type;


FieldMap::FieldMap(const std::string& fname_, size_t id_) :
		id(id_),
		nx(0), nz(0),
		x_min(0), dx(0), x_max(0),
		z_min(0), dz(0), z_max(0),
		data(0)
{
	bool consistent_dimensions;

	bool header = true;
	this->read_fieldmap_from_file(fname_, header, consistent_dimensions);

	if (!consistent_dimensions){
		header = false;
		this->read_fieldmap_from_file(fname_, header, consistent_dimensions);
	}

	// throws exception in case dimensions do not agree
	if (!consistent_dimensions) { throw FieldMapException::inconsistent_dimensions; }

	calc_interpolant();
}


void FieldMap::delete_data() {
	delete [] this->data;
}


void FieldMap::read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions)
{

	const size_t capacity_inc = 10 * 1000;

	// opens file
	std::ifstream file(fname_.c_str());
	if (!file.is_open()) throw FieldMapException::file_not_found;

	std::string word, line;
	if (header){
		// reads header section
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

		try{ this->physical_length = std::atof(comprimento.c_str())/1000.0; } catch (...) {}

	}

	std::getline(file, line);
	std::getline(file, line);

	// reads data section
	std::set<double> x_set, y_set, z_set;
	double x,y,z,bx,by,bz;
	size_t nr_points = 0;
	size_t capacity  = 0;
	if (this->data != 0) { std::free(this->data); } // deallocates memory if already set
	this->data = (double*) std::malloc(capacity*3*sizeof(double));

	while (true) {
		if (nr_points >= capacity) {
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
	this->nz     = z_set.size();
	this->x_min  = *(x_set.begin());
	this->x_max  = *(x_set.rbegin());
	this->z_min  = *(z_set.begin());
	this->z_max  = *(z_set.rbegin());
	this->y      = *(y_set.begin());
	this->dx     = (this->nx > 1) ? (this->x_max - this->x_min)/(this->nx - 1) : 0.0;
	this->dz     = (this->nz > 1) ? (this->z_max - this->z_min)/(this->nz - 1) : 0.0;
	this->x_grid.assign(x_set.begin(), x_set.end());
	this->z_grid.assign(z_set.begin(), z_set.end());

	// std::cout << "nr.points       : " << nr_points    << std::endl;
	// std::cout << "nr.points.x_set : " << this->nx     << std::endl;
	// std::cout << "nr.points.z_set : " << this->nz     << std::endl;
	// std::cout << "min.x           : " << this->x_min  << std::endl;
	// std::cout << "max.x           : " << this->x_max  << std::endl;
	// std::cout << "min.z           : " << this->z_min  << std::endl;
	// std::cout << "max.z           : " << this->z_max  << std::endl;
	// std::cout << "dx              : " << this->dx     << std::endl;
	// std::cout << "dz              : " << this->dz     << std::endl;

	if (y_set.size() != 1 || nr_points != (this->nx * this->nz )) {
		consistent_dimensions = false;
	} else { consistent_dimensions = true; }

	file.close();

}

void FieldMap::calc_interpolant(){

	int data_array_size = 3*this->nx*this->nz;
	alglib::real_1d_array x_array;
	alglib::real_1d_array z_array;
	alglib::real_1d_array data_array;
	x_array.setcontent(this->nx, &this->x_grid[0]);
	z_array.setcontent(this->nz, &this->z_grid[0]);
	data_array.setcontent(data_array_size, &this->data[0]);

	alglib::spline2dinterpolant interpolant;
	alglib::spline2dbuildbicubicv(x_array, this->nx, z_array, this->nz, data_array, 3, interpolant);
	this->interpolant = interpolant;
	this->delete_data();

}

Vector3D<double> FieldMap::field(const Vector3D<double>& pos) const{

	alglib::real_1d_array field;
	alglib::spline2dcalcv(this->interpolant, pos.x, pos.z, field);
	return Vector3D<double>(field[0], field[1], field[2]);

}
