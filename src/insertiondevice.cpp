#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <api.h>

InsertionDevice::InsertionDevice(FieldMapContainer& fieldmap_container){
  this->fieldmaps = fieldmap_container; this->type = 0; this->set_attributes();
}

InsertionDevice::InsertionDevice(FieldMap& fieldmap){
  this->fieldmaps = FieldMapContainer(fieldmap); this->type = 0; this->set_attributes();
}

InsertionDevice::InsertionDevice(std::vector<FieldMap>& fieldmaps){
  this->fieldmaps = FieldMapContainer(fieldmaps); this->type = 0; this->set_attributes();
}

InsertionDevice::InsertionDevice(CassetteContainer& cassette_container){
  this->cassettes = cassette_container; this->type = 1; this->set_attributes();
}

InsertionDevice::InsertionDevice(HalbachCassette& cassette){
  this->cassettes = CassetteContainer(cassette); this->type = 1; this->set_attributes();
}

InsertionDevice::InsertionDevice(std::vector<HalbachCassette>& cassettes){
  this->cassettes = CassetteContainer(cassettes); this->type = 1; this->set_attributes();
}

InsertionDevice::InsertionDevice(Vector3D<> (*function)(Vector3D<>)){
  this->field_function = function; this->type = 2;
}

InsertionDevice::InsertionDevice(const InsertionDevice &obj){
  this->fieldmaps = obj.fieldmaps;
  this->cassettes = obj.cassettes;
  this->field_function = obj.field_function;
  this->type = obj.type;
  this->x_min = obj.x_min;
  this->x_max = obj.x_max;
  this->y_min = obj.y_min;
  this->y_max = obj.y_max;
  this->z_min = obj.z_min;
  this->z_max = obj.z_max;
  this->physical_length = obj.physical_length;
}

void InsertionDevice::set_attributes(){
  if (this->type == 0){
    this->x_min = this->fieldmaps.x_min;
    this->x_max = this->fieldmaps.x_max;
    this->y_min = this->fieldmaps.y_min;
    this->y_max = this->fieldmaps.y_max;
    this->z_min = this->fieldmaps.z_min;
    this->z_max = this->fieldmaps.z_max;
    this->physical_length = this->fieldmaps.physical_length;
  } else if (this->type == 1){
    this->x_min = this->cassettes.get_x_min();
    this->x_max = this->cassettes.get_x_max();
    this->y_min = this->cassettes.get_y_min();
    this->y_max = this->cassettes.get_y_max();
    this->z_min = this->cassettes.get_z_min();
    this->z_max = this->cassettes.get_z_max();
    this->physical_length = this->cassettes.get_physical_length();
  }
}

Vector3D<double> InsertionDevice::field( Vector3D<double>& pos) {
  Vector3D<double> field;
  if (this->type == 0){
    field = this->fieldmaps.field(pos);
  } else if (this->type == 1){
    field = this->cassettes.field(pos);
  } else if (this->type == 2){
    field = this->field_function(pos);
  }
  return field;
}

std::vector<Vector3D<double> > InsertionDevice::field(std::vector<Vector3D<double> >& pos) {
  std::vector<Vector3D<> > field;
  for (int i=0; i < pos.size(); i+=1){
    field.push_back(this->field(pos[i]));
  }
  return field;
}

void InsertionDevice::calc_trajectory(double brho, double beta, double runge_kutta_step, Mask& mask, Vector3D<> r, Vector3D<> p, std::vector<std::vector<double> >& trajectory){
  runge_kutta(*this, brho, beta, runge_kutta_step, mask, r, p, trajectory);
}

void InsertionDevice::calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step, KickMap& kickmap){
  double beta; double brho;
  calc_brho(energy, brho, beta);

  int count = 0;
  int size = grid.nx * grid.ny;
  Vector3D<> r(0.0, 0.0, this->z_min);
  Vector3D<> p(0.0, 0.0, 1.0);
  Vector3D<> kick;

  std::cout << std::endl << "Calculating kickmap..." << std::endl;

  std::vector<double> kick_x_vector;
  std::vector<double> kick_y_vector;
  std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

  for(int i = 0; i < grid.ny; i+=1){
    for(int j =0; j < grid.nx; j+=1){
      r.x = grid.x[j]; r.y = grid.y[i];
      runge_kutta(*this, brho, beta, runge_kutta_step, mask, r, p, kick);
      kick_x_vector.push_back(kick.x );
      kick_y_vector.push_back(kick.y);

      count += 1;
      if (count%10 == 0) {
        std::cout << std::setw(6) << std::setprecision(4) << std::setfill(' ') << 100.0*(double(count)/double(size)) << '%' << '\r' << std::flush;
      }

    }
    kick_x.push_back(kick_x_vector);
    kick_y.push_back(kick_y_vector);
    kick_x_vector.clear();
    kick_y_vector.clear();
  }

  KickMap temp_kickmap(this->physical_length, grid.x, grid.y, kick_x, kick_y);
  kickmap = temp_kickmap;
}


void InsertionDevice::calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step, KickMap& kickmap, std::vector<std::vector<std::vector<double> > >& trajectories){
  double beta; double brho;
  calc_brho(energy, brho, beta);

  int count = 0;
  int size = grid.nx * grid.ny;
  Vector3D<> r(0.0, 0.0, this->z_min);
  Vector3D<> p(0.0, 0.0, 1.0);

  std::cout << std::endl << "Calculating kickmap..." << std::endl;

  std::vector<double> kick_x_vector;
  std::vector<double> kick_y_vector;
  std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;
  std::vector<std::vector<double> > trajectory;

  for(int i = 0; i < grid.ny; i+=1){
    for(int j =0; j < grid.nx; j+=1){
      r.x = grid.x[j]; r.y = grid.y[i];
      runge_kutta(*this, brho, beta, runge_kutta_step, mask, r, p, trajectory);
      kick_x_vector.push_back( (trajectory.back()[3])*(pow(brho, 2.0)) );
      kick_y_vector.push_back( (trajectory.back()[4])*(pow(brho, 2.0)) );
      trajectories.push_back(trajectory);

      count += 1;
      if (count%10 == 0) {
        std::cout << std::setw(6) << std::setprecision(4) << std::setfill(' ') << 100.0*(double(count)/double(size)) << '%' << '\r' << std::flush;
      }

    }
    kick_x.push_back(kick_x_vector);
    kick_y.push_back(kick_y_vector);
    kick_x_vector.clear();
    kick_y_vector.clear();
  }

  KickMap temp_kickmap(this->physical_length, grid.x, grid.y, kick_x, kick_y);
  kickmap = temp_kickmap;
}
