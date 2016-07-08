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
#include <api.h>

struct InputParameters{
  int nr_fieldmaps;
  std::vector<std::string> fieldmap_filenames;
  std::string kickmap_filename;
  double energy;
  double rkstep;
  int grid_nx;
  int grid_ny;
  double grid_xmin;
  double grid_xmax;
  double grid_ymin;
  double grid_ymax;
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
    if (nr_inputs == (13 + inputs.nr_fieldmaps)){
      inputs.mask_height  = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.mask_width   = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.mask_shape   =  values.back(); values.pop_back();
      status = true;
    } else if (nr_inputs == (11 + inputs.nr_fieldmaps)){
      inputs.mask_filename      = values.back(); values.pop_back();
      inputs.mask_shape_in_file = true;
      status = true;
    } else{
      std::cout << "Invalid number of input parameters!" << std::endl;
      status = false;
    }

    if (status){
      inputs.grid_ymax          = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_ymin          = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_xmax          = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_xmin          = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.grid_ny            = std::atoi(values.back().c_str());          values.pop_back();
      inputs.grid_nx            = std::atoi(values.back().c_str());          values.pop_back();
      inputs.rkstep             = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
      inputs.energy             = std::atof(values.back().c_str());          values.pop_back();
      inputs.kickmap_filename   = values.back();                             values.pop_back();
      for (int i=0; i< inputs.nr_fieldmaps; i+=1) {
        inputs.fieldmap_filenames.push_back(values.back()); values.pop_back();
      }
    }

    std::cout << std::endl;
    std::cout << "Inputs:" << std::endl;
    std::cout << "Number of fieldmaps:  " << inputs.nr_fieldmaps << std::endl;
    std::cout << "Fieldmap filenames:   " << std::endl;
    for (int i=0; i< inputs.nr_fieldmaps; i+=1) { std::cout << "  " << inputs.fieldmap_filenames[i] << std::endl; }
    std::cout << "Kickmap filename:     " << std::endl <<  "  " << inputs.kickmap_filename << std::endl;
    std::cout << "Energy:               " << inputs.energy << std::endl;
    std::cout << "Runge-kutta step      " << inputs.rkstep << std::endl;
    std::cout << "Grid nx:              " << inputs.grid_nx << std::endl;
    std::cout << "Grid ny:              " << inputs.grid_ny << std::endl;
    std::cout << "Grid xmin:            " << inputs.grid_xmin << std::endl;
    std::cout << "Grid xmax:            " << inputs.grid_xmax << std::endl;
    std::cout << "Grid ymin:            " << inputs.grid_ymin << std::endl;
    std::cout << "Grid ymax:            " << inputs.grid_ymax << std::endl;
    std::cout << "Mask shape:           " << inputs.mask_shape << std::endl;
    std::cout << "Mask width:           " << inputs.mask_width << std::endl;
    std::cout << "Mask height:          " << inputs.mask_height << std::endl;
    std::cout << "Mask filename:        " << inputs.mask_filename << std::endl;
    std::cout << std::endl;

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

    bool use_field_simmetry = false;
    FieldMapContainer fieldmaps(inputs.fieldmap_filenames, use_field_simmetry);
    InsertionDevice insertiondevice(fieldmaps);

    Grid grid(inputs.grid_nx, inputs.grid_ny, inputs.grid_xmin, inputs.grid_xmax, inputs.grid_ymin, inputs.grid_ymax);

    Mask mask;
    if (inputs.mask_shape_in_file){ mask.load(inputs.mask_filename); }
    else { mask.load(inputs.mask_shape, inputs.mask_width, inputs.mask_height); }

    KickMap kickmap(insertiondevice, grid, mask, inputs.energy, inputs.rkstep);
    kickmap.write_kickmap(inputs.kickmap_filename);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time: " << std::setprecision(4) << elapsed << " s" << std::endl;
    std::cout << std::endl;
  }

}
