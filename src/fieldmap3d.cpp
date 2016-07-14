#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <set>
#include <api.h>

FieldMap3D::FieldMap3D(const std::string& fname_, size_t id_) :
		id(id_),
		nx(0), nz(0), ny(0),
		x_min(0), dx(0), x_max(0),
		y_min(0), dy(0), y_max(0),
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
	if (!consistent_dimensions) { throw InsertionDeviceException::inconsistent_dimensions; }

	calc_interpolants();
}

FieldMap3D::FieldMap3D(const FieldMap3D &obj){
	this->fname  			 = obj.fname;
	this->data   			 = obj.data;
	this->nx     			 = obj.nx;
	this->ny     			 = obj.ny;
	this->nz     			 = obj.nz;
	this->x_min  			 = obj.x_min;
	this->x_max  			 = obj.x_max;
	this->y_min  			 = obj.y_min;
	this->y_max  			 = obj.y_max;
	this->z_min  			 = obj.z_min;
	this->z_max  			 = obj.z_max;
	this->dx     			 = obj.dx;
	this->dy      		 = obj.dy;
	this->dz     			 = obj.dz;
	this->x_grid 			 = obj.x_grid;
	this->y_grid 			 = obj.y_grid;
	this->z_grid 		   = obj.z_grid;
	this->interpolants = obj.interpolants;
}


void FieldMap3D::delete_data() {
	if (this->data != 0) { std::free(this->data); }
}

void FieldMap3D::read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions)
{

	const size_t capacity_inc = 10 * 1000;

	// opens file
	std::ifstream file(fname_.c_str());
	if (!file.is_open()) throw InsertionDeviceException::file_not_found;

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

	if (nr_points != (this->nx * this->nz * this->ny)) {
		consistent_dimensions = false;
	} else { consistent_dimensions = true; }

	file.close();

}

void FieldMap3D::get_y_data(std::vector<double>& y_data, int y_index){
	if (!y_data.empty()) { y_data.clear(); }

	int init = 3*y_index*this->nx;
	for (int j; j < this->nz; j+=1){
		for (int i; i <  this->nx; i+=1){
			y_data.push_back(this->data[init + 3*(i+j*this->nx*this->ny)+0]);
			y_data.push_back(this->data[init + 3*(i+j*this->nx*this->ny)+1]);
			y_data.push_back(this->data[init + 3*(i+j*this->nx*this->ny)+2]);
		}
	}

}

void FieldMap3D::calc_interpolants(){
	if (!this->interpolants.empty()) { this->interpolants.clear(); }

	std::vector<double> y_data;
	int data_y_array_size = 3*this->nx*this->nz;
	alglib::real_1d_array x_array;
	alglib::real_1d_array z_array;
	alglib::real_1d_array data_y_array;
	x_array.setcontent(this->nx, &this->x_grid[0]);
	z_array.setcontent(this->nz, &this->z_grid[0]);

	for (int i; i < this->ny; i+=1){
		this->get_y_data(y_data, i);
		data_y_array.setcontent(data_y_array_size, &y_data[0]);
		alglib::spline2dinterpolant interpolant;
		alglib::spline2dbuildbilinearv(x_array, this->nx, z_array, this->nz, data_y_array, 3, interpolant);
		this->interpolants.push_back(interpolant);
	}

}

Vector3D<double> FieldMap3D::field2D(const Vector3D<double>& pos, int interpolant_index) const{

	alglib::real_1d_array field;
	alglib::spline2dcalcv(this->interpolants[interpolant_index], pos.x, pos.z, field);
	return Vector3D<double>(field[0], field[1], field[2]);

}

Vector3D<double> FieldMap3D::field(const Vector3D<double>& pos) const{

	Vector3D<> field;

  if (this->ny == 1){

		field = this->field2D(pos, 0);

  } else {

    Vector3D<> f;
    std::vector<double> y_vector;
    std::vector<double> bx_vector;
    std::vector<double> by_vector;
    std::vector<double> bz_vector;

    for(int i=0; i < ny; i+=1){
      f = this->field2D(pos, i);
      bx_vector.push_back(f.x);
      by_vector.push_back(f.y);
      bz_vector.push_back(f.z);
    }

    alglib::real_1d_array y_array;
    alglib::real_1d_array bx_array;
    alglib::real_1d_array by_array;
    alglib::real_1d_array bz_array;
    y_array.setcontent(this->y_grid.size(), &this->y_grid[0]);
    bx_array.setcontent(bx_vector.size(), &bx_vector[0]);
    by_array.setcontent(by_vector.size(), &by_vector[0]);
    bz_array.setcontent(bz_vector.size(), &bz_vector[0]);

    alglib::spline1dinterpolant interpolant_x;
    alglib::spline1dinterpolant interpolant_y;
    alglib::spline1dinterpolant interpolant_z;
    alglib::spline1dbuildcubic(y_array, bx_array, interpolant_x);
    alglib::spline1dbuildcubic(y_array, by_array, interpolant_y);
    alglib::spline1dbuildcubic(y_array, bz_array, interpolant_z);

    field.x = alglib::spline1dcalc(interpolant_x, pos.y);
    field.y = alglib::spline1dcalc(interpolant_y, pos.y);
    field.z = alglib::spline1dcalc(interpolant_z, pos.y);
  }

  return field;

}

std::vector<Vector3D<double> > FieldMap3D::field(const std::vector<Vector3D<double> >& pos) const{

	std::vector<Vector3D<> > field;
  for (int i=0; i < pos.size(); i+=1){
    field.push_back(this->field(pos[i]));
  }
  return field;

}
