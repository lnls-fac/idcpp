#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <set>
#include <algorithm>
#include <fieldmap.h>
#include <stdlib.h>

static const double electron_rest_energy = 510998.92811;  // [eV]
static const double light_speed          = 299792458;     // [m/s]

struct InputParameters{
  std::vector<std::string> fieldmap_filenames;
  int nr_fieldmaps;
  std::string kickmap_filename;
  double energy;
  double rk_step;
  int nrpts_x;
  int nrpts_y;
  double width;
  double heigth;
};

void read_input_file(std::string input_filename, bool& status, InputParameters& inputs){
  std::ifstream input_file(input_filename.c_str());
  std::string line;
  std::string delimiter = " ";
  std::vector<std::string> values;

  if (!input_file.is_open()){
    std::cout << "Can't open input file!" << std::endl;
    status = false;
  }
  else{
    while(!input_file.eof()){
      getline(input_file, line);
      if (line.find("#")!=0 && !line.empty()) values.push_back(line);
    }
    inputs.heigth               = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
    inputs.width                = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
    inputs.nrpts_y              = std::atoi(values.back().c_str());          values.pop_back();
    inputs.nrpts_x              = std::atoi(values.back().c_str());          values.pop_back();
    inputs.rk_step              = std::atof(values.back().c_str());          values.pop_back();
    inputs.energy               = std::atof(values.back().c_str());          values.pop_back();
    inputs.kickmap_filename     = values.back();                             values.pop_back();
    inputs.nr_fieldmaps         = std::atoi(values.back().c_str());          values.pop_back();
    for (int i=0; i< inputs.nr_fieldmaps; i+=1) {
      inputs.fieldmap_filenames.push_back(values.back());
      values.pop_back();
    }
    status = true;
  }
}


void load_fieldmaps(InputParameters inputs, std::vector<FieldMap>& fieldmaps, bool& status){
  try{
    std::string filename;
    for(int i=0; i < inputs.nr_fieldmaps; i+=1){
      filename = inputs.fieldmap_filenames.back(); inputs.fieldmap_filenames.pop_back();
      FieldMap fieldmap(filename.c_str());
      fieldmaps.push_back(fieldmap);
      status = true;
      std::cout << "Load " << filename << std::endl;
    }
  } catch (...){
    status = false;
  }
}

void grid(int nrpts_x, int nrpts_y, double width, double heigth, std::vector<double>& x, std::vector<double>& y){
  if (nrpts_x > 1){
    for(int j = 0; j < nrpts_x; j+=1){
      x.push_back(-width + j*(2.0*width/(double(nrpts_x) - 1.0)));
    }
  } else x.push_back(0);
  if (nrpts_y > 1){
    for(int i = 0; i < nrpts_y; i+=1){
      y.push_back(- heigth + i*(2.0*heigth/(double(nrpts_y) - 1.0)));
    }
  } else y.push_back(0);
}

void calc_brho(double energy, double& beta, double& brho){
  double gamma = energy / electron_rest_energy;
  beta  = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  brho  = (beta * energy / light_speed);
}

void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds){
  dr_ds.x = p.x;
  dr_ds.y = p.y;
  dr_ds.z = p.z;
  dp_ds.x = - alpha * (p.y * b.z - p.z * b.y);
  dp_ds.y = - alpha * (p.z * b.x - p.x * b.z);
  dp_ds.z = - alpha * (p.x * b.y - p.y * b.x);
}

void get_interpolated_field3D(Vector3D<> r, std::vector<double> y, std::vector<double> bx, std::vector<double> by, std::vector<double> bz, Vector3D<>& b){
  alglib::real_1d_array y_array;
  y_array.setcontent(y.size(), &y[0]);

  alglib::real_1d_array bx_array;
  alglib::real_1d_array by_array;
  alglib::real_1d_array bz_array;
  bx_array.setcontent(bx.size(), &bx[0]);
  by_array.setcontent(by.size(), &by[0]);
  bz_array.setcontent(bz.size(), &bz[0]);

  alglib::spline1dinterpolant interpolant_x;
  alglib::spline1dinterpolant interpolant_y;
  alglib::spline1dinterpolant interpolant_z;
  alglib::spline1dbuildcubic(y_array, bx_array, interpolant_x);
  alglib::spline1dbuildcubic(y_array, by_array, interpolant_y);
  alglib::spline1dbuildcubic(y_array, bz_array, interpolant_z);

  b.x = alglib::spline1dcalc(interpolant_x, r.y);
  b.y = alglib::spline1dcalc(interpolant_y, r.y);
  b.z = alglib::spline1dcalc(interpolant_z, r.y);
}


void test_interpolation(InputParameters inputs){

  std::vector<FieldMap> fieldmaps; bool status;
  load_fieldmaps(inputs, fieldmaps, status);

  Vector3D<> r;
  Vector3D<> b;
  Vector3D<> b_calc;
  std::vector<double> diff;
  int nrpts = 50;

  Vector3D<> f;
  std::vector<double> bx_vector;
  std::vector<double> by_vector;
  std::vector<double> bz_vector;
  std::vector<double> y_values;
  for(int i=0; i < fieldmaps.size(); i+=1) y_values.push_back(fieldmaps[i].y);

  double xmax = fieldmaps[0].x_max;
  double xmin = fieldmaps[0].x_min;
  double ymin = *std::min_element(y_values.begin(), y_values.end());
  double ymax = *std::max_element(y_values.begin(), y_values.end());
  double zmax = fieldmaps[0].z_max;
  double zmin = fieldmaps[0].z_min;
  double dx = (xmax - xmin)/(double(nrpts)-1.0);
  double dy = (ymax - ymin)/(double(nrpts)-1.0);
  double dz = (zmax - zmin)/(double(nrpts)-1.0);

  for (int i=0; i < nrpts; i+=1){
    for (int j=0; j < nrpts; j+=1){
      for (int k=0; k < nrpts; k+=1){
        r.x = xmin + i*dx;
        r.y = ymin + j*dy;
        r.z = zmin + k*dz;

        for(int i=0; i < fieldmaps.size(); i+=1){
          try{ f = fieldmaps[i].field(r); }
          catch (...) { f.x = f.y = f.z = 0.0; }
          bx_vector.push_back(f.x);
          by_vector.push_back(f.y);
          bz_vector.push_back(f.z);
        }
        get_interpolated_field3D(r, y_values, bx_vector, by_vector, bz_vector, b);
        bx_vector.clear(); by_vector.clear(); bz_vector.clear();

        b_calc.x = r.z*r.z;
        b_calc.y = r.x*r.x + r.y;
        b_calc.z = r.z*r.y*r.x;
        diff.push_back((b.x - b_calc.x));
        diff.push_back((b.y - b_calc.y));
        diff.push_back((b.z - b_calc.z));

        std::cout <<  r << " " <<(b.x - b_calc.x) << " " << (b.y - b_calc.y) << " " << (b.z - b_calc.z) << std::endl;
      }
    }
  }

  std::cout << std::endl;
  std::cout << "Test interpolation - Max diff: "<<*std::max_element(diff.begin(), diff.end()) << std::endl;
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

  if(status) {
    test_interpolation(inputs);
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time: " << elapsed << " s" << std::endl;
    std::cout << std::endl;
  }

}
