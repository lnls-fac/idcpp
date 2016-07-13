#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <api.h>

KickMap::KickMap(double physical_length, std::vector<double> x, std::vector<double> y, std::vector<std::vector<double> > kick_x, std::vector<std::vector<double> > kick_y){
  this->physical_length = physical_length;
  this->nx = x.size();
  this->ny = y.size();
  this->x = x;
  this->y = y;
  this->kick_x = kick_x;
  this->kick_y = kick_y;
}

void KickMap::write_to_file(std::string filename){
  _write_kickmap_file(filename, this->physical_length, this->x, this->y, this->kick_x, this->kick_y);
}


InsertionDevice::InsertionDevice(FieldMapContainer& fieldmap_container){
  this->fieldmaps = fieldmap_container;
  this->x_min = this->fieldmaps.x_min;
  this->x_max = this->fieldmaps.x_max;
  this->y_min = this->fieldmaps.y_min;
  this->y_max = this->fieldmaps.y_max;
  this->z_min = float(this->fieldmaps.z_min);
  this->z_max = float(this->fieldmaps.z_max);
  this->physical_length = this->fieldmaps.physical_length;
  this->type = 0;
}

InsertionDevice::InsertionDevice(FieldMap& fieldmap){
  FieldMapContainer fieldmap_container(fieldmap);
  this->fieldmaps = fieldmap_container;
  this->x_min = this->fieldmaps.x_min;
  this->x_max = this->fieldmaps.x_max;
  this->y_min = this->fieldmaps.y_min;
  this->y_max = this->fieldmaps.y_max;
  this->z_min = this->fieldmaps.z_min;
  this->z_max = this->fieldmaps.z_max;
  this->physical_length = this->fieldmaps.physical_length;
  this->type = 0;
}

InsertionDevice::InsertionDevice(std::vector<FieldMap>& fieldmaps){
  FieldMapContainer fieldmap_container(fieldmaps);
  this->fieldmaps = fieldmap_container;
  this->x_min = this->fieldmaps.x_min;
  this->x_max = this->fieldmaps.x_max;
  this->y_min = this->fieldmaps.y_min;
  this->y_max = this->fieldmaps.y_max;
  this->z_min = this->fieldmaps.z_min;
  this->z_max = this->fieldmaps.z_max;
  this->physical_length = this->fieldmaps.physical_length;
  this->type = 0;
}

InsertionDevice::InsertionDevice(EPU& epu){
  this->epu = epu;
  this->x_min = this->epu.get_x_min();
  this->x_max = this->epu.get_x_max();
  this->y_min = this->epu.get_y_min();
  this->y_max = this->epu.get_y_max();
  this->z_min = this->epu.get_z_min();
  this->z_max = this->epu.get_z_max();
  this->physical_length = this->epu.get_physical_length();
  this->type = 1;
}

InsertionDevice::InsertionDevice(DELTA& delta){
  this->delta = delta;
  this->x_min = this->delta.get_x_min();
  this->x_max = this->delta.get_x_max();
  this->y_min = this->delta.get_y_min();
  this->y_max = this->delta.get_y_max();
  this->z_min = this->delta.get_z_min();
  this->z_max = this->delta.get_z_max();
  this->physical_length = this->delta.get_physical_length();
  this->type = 2;
}

InsertionDevice::InsertionDevice(Vector3D<> (*function)(Vector3D<>)){
  this->field_function = function;
  this->type = 3;
}

Vector3D<double> InsertionDevice::field( Vector3D<double>& pos) {
  Vector3D<double> field;
  if (this->type == 0){
    field = this->fieldmaps.field(pos);
  } else if (this->type == 1){
    field = this->epu.field(pos);
  } else if (this->type == 2){
    field = this->delta.field(pos);
  } else if (this->type == 3){
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

void InsertionDevice::calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step){
  KickMap kickmap;
  std::vector<std::vector<std::vector<double> > > trajectories;
  _calc_kickmap(*this, grid, mask, energy, runge_kutta_step, this->z_min, this->z_max, this->physical_length, kickmap, trajectories);
  this->kickmap = kickmap;
}

void InsertionDevice::calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step, std::string trajectories_filename){
  KickMap kickmap;
  std::vector<std::vector<std::vector<double> > > trajectories;
  _calc_kickmap(*this, grid, mask, energy, runge_kutta_step, this->z_min, this->z_max, this->physical_length, kickmap, trajectories);
  this->kickmap = kickmap;
  save_trajectories(trajectories_filename, trajectories);
}

void InsertionDevice::write_kickmap_file(std::string filename){
  this->kickmap.write_to_file(filename);
}

void InsertionDevice::write_fieldmap_files(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  _write_fieldmap_files(*this, filename, x_vector, y_vector, z_vector);
}

void InsertionDevice::write_fieldmap_file(std::string filename, std::vector<double> x_vector, double y, std::vector<double> z_vector){
  std::vector<double> y_vector;
  y_vector.push_back(y);
  _write_fieldmap_file(*this, filename, x_vector, y_vector, z_vector);
}

void InsertionDevice::write_fieldmap_file(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  _write_fieldmap_file(*this, filename, x_vector, y_vector, z_vector);
}
