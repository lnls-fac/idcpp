#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
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
  output_file << this->physical_length << std::endl;
  output_file << "# Number of Horizontal Points" << std::endl;
  output_file << this->nx << std::endl;
  output_file << "# Number of Vertical Points" << std::endl;
  output_file << this->ny << std::endl;
  output_file << "# Horizontal KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";

  for(int j =0; j < this->nx; j+=1){
      output_file << std::scientific << std::showpos << this->x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < this->ny; i+=1){
    output_file << std::scientific << std::showpos << this->y[i] << " ";
    for(int j =0; j < this->nx; j+=1){
      if ( isnan(this->kick_x[i][j]) ) { output_file << std::setw(13) << std::setfill(' ') << this->kick_x[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << this->kick_x[i][j] << " "; }
    }
    output_file << std::endl;
  }

  output_file << "# Vertical KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";
  for(int j =0; j < this->nx; j+=1){
      output_file << std::scientific << std::showpos << this->x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < this->ny; i+=1){
    output_file << std::scientific << std::showpos << this->y[i] << " ";
    for(int j =0; j < this->nx; j+=1){
      if ( isnan(this->kick_y[i][j]) ){ output_file << std::setw(13) << std::setfill(' ') << this->kick_y[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << this->kick_y[i][j] << " "; }
    }
    output_file << std::endl;
  }

  std::cout << "Kickmap saved in file: " << std::endl;
  std::cout << filename << std::endl;
  std::cout << std::endl;

}
