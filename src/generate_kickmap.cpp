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
  double height;
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
      getline(input_file, line);
      if (line.find("#")!=0 && !line.empty()) values.push_back(line);
    }
    inputs.height               = (std::atof(values.back().c_str()))/1000.0; values.pop_back();
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
  std::cout << std::endl;
  std::cout << "Loading fieldmap files..." << std::endl;
  std::string filename;
  try{
    for(int i=0; i < inputs.nr_fieldmaps; i+=1){
      filename = inputs.fieldmap_filenames.back(); inputs.fieldmap_filenames.pop_back();
      FieldMap fieldmap(filename.c_str());
      fieldmaps.push_back(fieldmap);
      status = true;
      std::cout << "Load " << filename << std::endl;
    }
  } catch (...){
    std::cout << "Can't open fieldmap file: " << filename << std::endl;
    std::cout << std::endl;
    status = false;
  }
}

void grid(int nrpts_x, int nrpts_y, double width, double height, std::vector<double>& x, std::vector<double>& y){
  if (nrpts_x > 1){
    for(int j = 0; j < nrpts_x; j+=1){
      x.push_back(-width + j*(2.0*width/(double(nrpts_x) - 1.0)));
    }
  } else x.push_back(0);
  if (nrpts_y > 1){
    for(int i = 0; i < nrpts_y; i+=1){
      y.push_back(- height + i*(2.0*height/(double(nrpts_y) - 1.0)));
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

void get_interpolated_field3D(Vector3D<> r, std::vector<double> y, std::vector<double>& bx, std::vector<double>& by, std::vector<double>& bz, Vector3D<>& b){
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

  bx.clear(); by.clear(); bz.clear();
}

void check_position(InputParameters inputs, Vector3D<> r, bool& inside){
  double ymax;
  double x = abs(r.x);
  double y = abs(r.y);
  ymax = inputs.width + (inputs.width/inputs.height)*r.x;
  inside = (y < ymax) ? true : false;
}

void runge_kutta(InputParameters inputs, std::vector<FieldMap> fieldmaps, double brho, double beta, double max_rz, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks){

  double alpha = 1.0/brho/beta;
  double s_step = inputs.rk_step;
  bool inside = true;
  Vector3D<> b; Vector3D<> b1; Vector3D<> b2; Vector3D<> b3;
  Vector3D<> kr1; Vector3D<> kp1; Vector3D<> r1; Vector3D<> p1;
  Vector3D<> kr2; Vector3D<> kp2; Vector3D<> r2; Vector3D<> p2;
  Vector3D<> kr3; Vector3D<> kp3; Vector3D<> r3; Vector3D<> p3;
  Vector3D<> kr4; Vector3D<> kp4;

  if (fieldmaps.size() == 1){

    while (r.z < max_rz){

      check_position(inputs, r, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      b = fieldmaps[0].field(r);
      newton_lorentz_equation(alpha, r, p, b, kr1, kp1);
      r1 = r + (s_step/2.0)* kr1;
      p1 = p + (s_step/2.0)* kp1;

      check_position(inputs, r1, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      b1 = fieldmaps[0].field(r1);
      newton_lorentz_equation(alpha, r1, p1, b1, kr2, kp2);
      r2 = r + (s_step/2.0)* kr2;
      p2 = p + (s_step/2.0)* kp2;

      check_position(inputs, r2, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      b2 = fieldmaps[0].field(r2);
      newton_lorentz_equation(alpha, r2, p2, b2, kr3, kp3);
      r3 = r + s_step* kr3;
      p3 = p + s_step* kp3;

      check_position(inputs, r3, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      b3 = fieldmaps[0].field(r3);
      newton_lorentz_equation(alpha, r3, p3, b3, kr4, kp4);

      r = r + (s_step/6.0)*(kr1 + 2.0*kr2 + 2.0*kr3 + kr4);
      p = p + (s_step/6.0)*(kp1 + 2.0*kp2 + 2.0*kp3 + kp4);

    }

  } else {

    Vector3D<> f;
    std::vector<double> bx_vector;
    std::vector<double> by_vector;
    std::vector<double> bz_vector;
    std::vector<double> y_values;
    for(int i=0; i < fieldmaps.size(); i+=1) y_values.push_back(fieldmaps[i].y);

    while (r.z < max_rz){

      check_position(inputs, r, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      for(int i=0; i < fieldmaps.size(); i+=1){
        f = fieldmaps[i].field(r);
        bx_vector.push_back(f.x); by_vector.push_back(f.y); bz_vector.push_back(f.z);
      }
      get_interpolated_field3D(r, y_values, bx_vector, by_vector, bz_vector, b);
      newton_lorentz_equation(alpha, r, p, b, kr1, kp1);
      r1 = r + (s_step/2.0)* kr1;
      p1 = p + (s_step/2.0)* kp1;

      check_position(inputs, r1, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      for(int i=0; i < fieldmaps.size(); i+=1){
        f = fieldmaps[i].field(r1);
        bx_vector.push_back(f.x); by_vector.push_back(f.y); bz_vector.push_back(f.z);
      }
      get_interpolated_field3D(r1, y_values, bx_vector, by_vector, bz_vector, b1);
      newton_lorentz_equation(alpha, r1, p1, b1, kr2, kp2);
      r2 = r + (s_step/2.0)* kr2;
      p2 = p + (s_step/2.0)* kp2;

      check_position(inputs, r2, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      for(int i=0; i < fieldmaps.size(); i+=1){
        f = fieldmaps[i].field(r2);
        bx_vector.push_back(f.x); by_vector.push_back(f.y); bz_vector.push_back(f.z);
      }
      get_interpolated_field3D(r2, y_values, bx_vector, by_vector, bz_vector, b2);
      newton_lorentz_equation(alpha, r2, p2, b2, kr3, kp3);
      r3 = r + s_step* kr3;
      p3 = p + s_step* kp3;

      check_position(inputs, r3, inside); if(!inside) { p.x = p.y = p.z = NAN; break; }
      for(int i=0; i < fieldmaps.size(); i+=1){
        f = fieldmaps[i].field(r3);
        bx_vector.push_back(f.x); by_vector.push_back(f.y); bz_vector.push_back(f.z);
      }
      get_interpolated_field3D(r3, y_values, bx_vector, by_vector, bz_vector, b3);
      newton_lorentz_equation(alpha, r3, p3, b3, kr4, kp4);

      r = r + (s_step/6.0)*(kr1 + 2.0*kr2 + 2.0*kr3 + kr4);
      p = p + (s_step/6.0)*(kp1 + 2.0*kp2 + 2.0*kp3 + kp4);

    }
  }
  kicks = p*(pow(brho, 2.0));
}

void generate_kickmap(InputParameters inputs, bool& status){

  time_t now = time(0); char* date = ctime(&now);

  std::vector<FieldMap> fieldmaps;
  load_fieldmaps(inputs, fieldmaps, status);

  if (status){
    std::ofstream output_file(inputs.kickmap_filename.c_str());

    if(! output_file) {
      status = false;
      std::cout << std::endl;
      std::cout << "Can't open output file: " << inputs.kickmap_filename << std::endl;
      std::cout << std::endl;
    } else {

      std::vector<double> x_grid; std::vector<double> y_grid;
      grid(inputs.nrpts_x, inputs.nrpts_y, inputs.width, inputs.height, x_grid, y_grid);

      double beta; double brho;
      calc_brho(inputs.energy, beta, brho);

      std::vector<double> z_min_vector;
      std::vector<double> z_max_vector;
      for(int i=0; i < fieldmaps.size(); i+=1) {
        z_min_vector.push_back(fieldmaps[i].z_min);
        z_max_vector.push_back(fieldmaps[i].z_max);
      }
      double min_rz = *std::min_element(z_min_vector.begin(), z_min_vector.end());
      double max_rz = *std::max_element(z_max_vector.begin(), z_max_vector.end());

      double kick_x[y_grid.size()][x_grid.size()];
      double kick_y[y_grid.size()][x_grid.size()];
      Vector3D<> r(0.0, 0.0, min_rz);
      Vector3D<> p(0.0, 0.0, 1.0);
      Vector3D<> kicks;

      std::cout << std::endl;
      std::cout << "Calculating kickmap..." << std::endl;
      int count = 0;
      int size = x_grid.size()*y_grid.size();
      for(int i = 0; i < y_grid.size(); i+=1){
        for(int j =0; j < x_grid.size(); j+=1){
          r.x = x_grid[j];
          r.y = y_grid[i];
          runge_kutta(inputs, fieldmaps, brho, beta, max_rz, r, p, kicks);
          kick_x[i][j] = kicks.x;
          kick_y[i][j] = kicks.y;
          count += 1;
          if (count%10 == 0) {
            std::cout << std::setw(6) << std::setprecision(4) << std::setfill(' ') << 100.0*(double(count)/double(size)) << '%' << '\r' << std::flush;
          }
        }
      }

      output_file << "# KICKMAP" << std::endl;
      output_file << "# Author: Luana N. P. Vilela @ LNLS, Date: " << date;
      output_file << "# ID Length [m]" << std::endl;
      output_file << fieldmaps[0].physical_length << std::endl;
      output_file << "# Number of Horizontal Points" << std::endl;
      output_file << inputs.nrpts_x << std::endl;
      output_file << "# Number of Vertical Points" << std::endl;
      output_file << inputs.nrpts_y << std::endl;
      output_file << "# Horizontal KickTable in T2m2" << std::endl;
      output_file << "START" << std::endl;
      output_file << std::setw(14) << " ";

      for(int j =0; j < x_grid.size(); j+=1){
          output_file << std::scientific << std::showpos << x_grid[j] << " ";
      }
      output_file << std::endl;
      for(int i = 0; i < y_grid.size(); i+=1){
        output_file << std::scientific << std::showpos << y_grid[i] << " ";
        for(int j =0; j < x_grid.size(); j+=1){
          output_file << std::scientific << std::showpos << kick_x[i][j] << " ";
        }
        output_file << std::endl;
      }

      output_file << "# Vertical KickTable in T2m2" << std::endl;
      output_file << "START" << std::endl;
      output_file << std::setw(14) << " ";
      for(int j =0; j < x_grid.size(); j+=1){
          output_file << std::scientific << std::showpos << x_grid[j] << " ";
      }
      output_file << std::endl;
      for(int i = 0; i < y_grid.size(); i+=1){
        output_file << std::scientific << std::showpos << y_grid[i] << " ";
        for(int j =0; j < x_grid.size(); j+=1){
          output_file << std::scientific << std::showpos << kick_y[i][j] << " ";
        }
        output_file << std::endl;
      }
      std::cout << "Kickmap saved in file: " << std::endl;
      std::cout << inputs.kickmap_filename << std::endl;
      std::cout << std::endl;
    }
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

  if(status) generate_kickmap(inputs, status);

  if (status){
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time: " << std::setprecision(4) << elapsed << " s" << std::endl;
    std::cout << std::endl;
  }

}
