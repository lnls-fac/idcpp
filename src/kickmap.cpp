#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <api.h>

static const double electron_rest_energy = 510998.92811;  // [eV]
static const double light_speed          = 299792458;     // [m/s]

KickMap::KickMap(std::vector<FieldMap> fieldmaps, Grid grid, Mask mask, double energy, double runge_kutta_step){
  this->fieldmaps = fieldmaps;
  this->grid = grid;
  this->mask = mask;
  this->energy = energy;
  this->runge_kutta_step = runge_kutta_step;
  this->calc_kicks();
}

void calc_brho(double energy, double& brho, double& beta){
  double gamma = energy / electron_rest_energy;
  beta  = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  brho  = (beta * energy / light_speed);
}

double get_zmin(std::vector<FieldMap> fieldmaps){
  std::vector<double> zmin_vector;
  for(int i=0; i < fieldmaps.size(); i+=1) {
    zmin_vector.push_back(fieldmaps[i].z_min);
  }
  double zmin = *std::min_element(zmin_vector.begin(), zmin_vector.end());
  return zmin;
}

double get_zmax(std::vector<FieldMap> fieldmaps){
  std::vector<double> zmax_vector;
  for(int i=0; i < fieldmaps.size(); i+=1) {
    zmax_vector.push_back(fieldmaps[i].z_max);
  }
  double zmax = *std::min_element(zmax_vector.begin(), zmax_vector.end());
  return zmax;
}

void KickMap::calc_kicks(){

  double beta; double brho;
  calc_brho(this->energy, brho, beta);

  double zmin = get_zmin(this->fieldmaps);
  double zmax = get_zmax(this->fieldmaps);

  int count = 0;
  int size = this->grid.nx * this->grid.ny;
  std::vector<double> kick_x_vector;
  std::vector<double> kick_y_vector;
  Vector3D<> r(0.0, 0.0, zmin);
  Vector3D<> p(0.0, 0.0, 1.0);
  Vector3D<> kick;

  std::cout << std::endl;
  std::cout << "Calculating kickmap..." << std::endl;

  this->kick_x.clear();
  this->kick_y.clear();

  for(int i = 0; i < this->grid.ny; i+=1){
    kick_x_vector.clear();
    kick_y_vector.clear();

    for(int j =0; j < this->grid.nx; j+=1){
      r.x = this->grid.x[j];
      r.y = this->grid.y[i];
      runge_kutta(this->fieldmaps, brho, beta, zmax, this->runge_kutta_step, this->mask, r, p, kick);
      kick_x_vector.push_back(kick.x);
      kick_y_vector.push_back(kick.y);

      count += 1;
      if (count%10 == 0) {
        std::cout << std::setw(6) << std::setprecision(4) << std::setfill(' ') << 100.0*(double(count)/double(size)) << '%' << '\r' << std::flush;
      }

    }
    this->kick_x.push_back(kick_x_vector);
    this->kick_y.push_back(kick_y_vector);
  }

}

void KickMap::write_kickmap(std::string filename){

  time_t now = time(0); char* date = ctime(&now);

  std::ofstream output_file(filename.c_str());
  if(!output_file) {
    std::cout << std::endl;
    std::cout << "Can't open output file: " << filename << std::endl;
    std::cout << std::endl;
    return;
  }

  output_file << "# KICKMAP" << std::endl;
  output_file << "# Author: Luana N. P. Vilela @ LNLS, Date: " << date;
  output_file << "# ID Length [m]" << std::endl;
  output_file << this->fieldmaps[0].physical_length << std::endl;
  output_file << "# Number of Horizontal Points" << std::endl;
  output_file << this->grid.nx << std::endl;
  output_file << "# Number of Vertical Points" << std::endl;
  output_file << this->grid.ny << std::endl;
  output_file << "# Horizontal KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";

  for(int j =0; j < this->grid.nx; j+=1){
      output_file << std::scientific << std::showpos << this->grid.x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < this->grid.ny; i+=1){
    output_file << std::scientific << std::showpos << this->grid.y[i] << " ";
    for(int j =0; j < this->grid.nx; j+=1){
      if ( isnan(this->kick_x[i][j]) ) { output_file << std::setw(13) << std::setfill(' ') << this->kick_x[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << this->kick_x[i][j] << " "; }
    }
    output_file << std::endl;
  }

  output_file << "# Vertical KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";
  for(int j =0; j < this->grid.nx; j+=1){
      output_file << std::scientific << std::showpos << this->grid.x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < this->grid.ny; i+=1){
    output_file << std::scientific << std::showpos << this->grid.y[i] << " ";
    for(int j =0; j < this->grid.nx; j+=1){
      if ( isnan(this->kick_y[i][j]) ){ output_file << std::setw(13) << std::setfill(' ') << this->kick_y[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << this->kick_y[i][j] << " "; }
    }
    output_file << std::endl;
  }

  std::cout << "Kickmap saved in file: " << std::endl;
  std::cout << filename << std::endl;
  std::cout << std::endl;

}
