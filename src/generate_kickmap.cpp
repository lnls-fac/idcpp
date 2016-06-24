#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <set>
#include <algorithm>
#include <stdlib.h>
#include <fieldmap.h>
#include <api.h>

struct InputParameters{
  int nr_fieldmaps;
  std::vector<std::string> fieldmap_filenames;
  std::string kickmap_filename;
  double energy;
  double rkstep;
  int grid_nx;
  int grid_ny;
  double grid_xwindow;
  double grid_ywindow;
  std::string mask_shape;
  double mask_width;
  double mask_height;
  std::string mask_filename;
  bool mask_shape_in_file = false;
};

void read_input_file(std::string input_filename, bool& status, InputParameters& inputs){
  std::ifstream input_file(input_filename.c_str());
  std::string line;
  std::string delimiter = " ";
  std::vector<std::string> values;

  if (!input_file.is_open()){
    std::cout << "Can't open input file: " << input_filename << std::endl;
    status = false;
  }
  else{
    while(!input_file.eof()){
      getline(input_file, line); if (line.find("#")!=0 && !line.empty()) { values.push_back(line);}
    }
    int nr_inputs = values.size();

    inputs.nr_fieldmaps = std::atoi(values.front().c_str());
    if (nr_inputs == (11 + inputs.nr_fieldmaps)){
      inputs.mask_height  = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.mask_width   = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.mask_shape   =  values.back(); values.pop_back();
      status = true;
    } else if (nr_inputs == (9 + inputs.nr_fieldmaps)){
      inputs.mask_filename      = values.back(); values.pop_back();
      inputs.mask_shape_in_file = true;
      status = true;
    } else{
      std::cout << "Invalid number of input parameters!" << std::endl;
      status = false;
    }

    if (status){
      inputs.grid_ywindow     = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_xwindow     = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_ny          = std::atoi(values.back().c_str());          values.pop_back();
      inputs.grid_nx          = std::atoi(values.back().c_str());          values.pop_back();
      inputs.rkstep           = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.energy           = std::atof(values.back().c_str());          values.pop_back();
      inputs.kickmap_filename = values.back();                             values.pop_back();
      for (int i=0; i< inputs.nr_fieldmaps; i+=1) {
        inputs.fieldmap_filenames.push_back(values.back()); values.pop_back();
      }
    }

  }
}

void load_fieldmaps(std::vector<std::string> fieldmap_filenames, std::vector<FieldMap>& fieldmaps, bool& status){
  std::cout << std::endl;
  std::cout << "Loading fieldmap files..." << std::endl;
  std::string filename;

  int nr_fieldmaps = fieldmap_filenames.size();
  try{
    for(int i=0; i < nr_fieldmaps; i+=1){
      filename = fieldmap_filenames.back();
      fieldmap_filenames.pop_back();
      FieldMap fieldmap(filename.c_str());
      fieldmaps.push_back(fieldmap);
      std::cout << "Load " << filename << std::endl;
      status = true;
    }
  } catch (...){
    std::cout << "Can't open fieldmap file: " << filename << std::endl;
    std::cout << std::endl;
    status = false;
  }
}

int main(int argc, char ** argv) {

  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);

  InputParameters inputs;
  bool status;

  if(argc != 2){
    std::cout << "Invalid number of arguments!" << std::endl;
    status = false;
  }
  else{
    std::string input_filename = argv[1];
    read_input_file(input_filename, status, inputs);
  }

  if (status){

    std::vector<FieldMap> fieldmaps;
    load_fieldmaps(inputs.fieldmap_filenames, fieldmaps, status);

    Grid grid(inputs.grid_nx, inputs.grid_ny, (inputs.grid_xwindow/2.0), (inputs.grid_ywindow/2.0));

    Mask mask;
    if (inputs.mask_shape_in_file){ mask.load(inputs.mask_filename); }
    else { mask.load(inputs.mask_shape, inputs.mask_width, inputs.mask_height); }

    KickMap kickmap(fieldmaps, grid, mask, inputs.energy, inputs.rkstep);
    kickmap.write_kickmap(inputs.kickmap_filename);

  }

  if (status){
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time: " << std::setprecision(4) << elapsed << " s" << std::endl;
    std::cout << std::endl;
  }

}
