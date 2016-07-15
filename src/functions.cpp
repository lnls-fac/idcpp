#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <api.h>

static const double electron_rest_energy = 510998.92811;  // [eV]
static const double light_speed          = 299792458;     // [m/s]


void create_epu(Block& genblock, unsigned int nr_periods, double magnetic_gap, double cassette_separation, double block_separation, double phase_csd, double phase_cie, InsertionDevice& insertiondevice){

  HalbachCassette csd; HalbachCassette cse; HalbachCassette cid; HalbachCassette cie;
  std::vector<HalbachCassette> cassettes;
  Vector3D<> mag = genblock.get_mag();
  if ((mag.x ==0) && (mag.y == 0) && (mag.z != 0)){
    csd.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
    cse.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
    genblock.set_mag(Vector3D<>(mag.x, mag.y, -mag.z));
    cid.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
    cie.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  } else if ((mag.x ==0) && (mag.y != 0) && (mag.z == 0)){
    csd.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
    cse.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90p(), nr_periods, block_separation);
    cid.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
    cie.gen_halbach_cassette(genblock,  Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  } else {
    throw InsertionDeviceException::invalid_magnetization;
  }

  Vector3D<double> dim = genblock.get_dim();
  csd.set_xcenter(+(cassette_separation + dim.x)/2.0);
  csd.set_ycenter(+(magnetic_gap + dim.y)/2.0);
  csd.set_zcenter(phase_csd);

  cse.set_xcenter(-(cassette_separation + dim.x)/2.0);
  cse.set_ycenter(+(magnetic_gap + dim.y)/2.0);
  cse.set_zcenter(0.0);

  cid.set_xcenter(+(cassette_separation + dim.x)/2.0);
  cid.set_ycenter(-(magnetic_gap + dim.y)/2.0);
  cid.set_zcenter(0.0);

  cie.set_xcenter(-(cassette_separation + dim.x)/2.0);
  cie.set_ycenter(-(magnetic_gap + dim.y)/2.0);
  cie.set_zcenter(phase_cie);

  cassettes.push_back(csd); cassettes.push_back(cse); cassettes.push_back(cid); cassettes.push_back(cie);
  insertiondevice = InsertionDevice(cassettes);
}


void create_delta(Block& genblock, unsigned int nr_periods, double vertical_gap, double horizontal_gap, double block_separation, double phase_cd, double phase_ce, InsertionDevice& insertiondevice){

  HalbachCassette cs; HalbachCassette cd; HalbachCassette ci; HalbachCassette ce;
  std::vector<HalbachCassette> cassettes;
  cs.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  cd.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90n(), nr_periods, block_separation);
  ci.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90p(), nr_periods, block_separation);
  ce.gen_halbach_cassette(genblock, Matrix3D<double>::rotx90p(), nr_periods, block_separation);

  Vector3D<double> dim = genblock.get_dim();
  cs.set_xcenter(0.0);
  cs.set_ycenter(+(vertical_gap + dim.y)/2.0);
  cs.set_zcenter(0.0);

  cd.set_xcenter(+(horizontal_gap + dim.x)/2.0);
  cd.set_ycenter(0.0);
  cd.set_zcenter(phase_cd);

  ci.set_xcenter(0.0);
  ci.set_ycenter(-(vertical_gap + dim.y)/2.0);
  ci.set_zcenter(0.0);

  ce.set_xcenter(-(horizontal_gap + dim.x)/2.0);
  ce.set_ycenter(0.0);
  ce.set_zcenter(phase_ce);

  cassettes.push_back(cs); cassettes.push_back(cd); cassettes.push_back(ci); cassettes.push_back(ce);
  insertiondevice = InsertionDevice(cassettes);
}


void calc_brho(double energy, double& brho, double& beta){
  double gamma = energy / electron_rest_energy;
  beta  = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  brho  = (beta * energy / light_speed);
}

void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds){

  dr_ds.x = p.x;
  dr_ds.y = p.y;
  dr_ds.z = p.z;
  dp_ds.x = -alpha * (p.y * b.z - p.z * b.y);
  dp_ds.y = -alpha * (p.z * b.x - p.x * b.z);
  dp_ds.z = -alpha * (p.x * b.y - p.y * b.x);

}


void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double step, Mask& mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kick){

  double alpha = 1.0/brho/beta;
  bool inside = true;
  Vector3D<> b; Vector3D<> b1; Vector3D<> b2; Vector3D<> b3;
  Vector3D<> kr1; Vector3D<> kp1; Vector3D<> r1; Vector3D<> p1;
  Vector3D<> kr2; Vector3D<> kp2; Vector3D<> r2; Vector3D<> p2;
  Vector3D<> kr3; Vector3D<> kp3; Vector3D<> r3; Vector3D<> p3;
  Vector3D<> kr4; Vector3D<> kp4;

  while (r.z < insertiondevice.z_max){

    inside = mask.is_inside(r); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b = insertiondevice.field(r);
    newton_lorentz_equation(alpha, r, p, b, kr1, kp1);
    r1 = r + (step/2.0)* kr1;
    p1 = p + (step/2.0)* kp1;

    inside = mask.is_inside(r1); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b1 = insertiondevice.field(r1);
    newton_lorentz_equation(alpha, r1, p1, b1, kr2, kp2);
    r2 = r + (step/2.0)* kr2;
    p2 = p + (step/2.0)* kp2;

    inside = mask.is_inside(r2); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b2 = insertiondevice.field(r2);
    newton_lorentz_equation(alpha, r2, p2, b2, kr3, kp3);
    r3 = r + step* kr3;
    p3 = p + step* kp3;

    inside = mask.is_inside(r3); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b3 = insertiondevice.field(r3);
    newton_lorentz_equation(alpha, r3, p3, b3, kr4, kp4);

    r = r + (step/6.0)*(kr1 + 2.0*kr2 + 2.0*kr3 + kr4);
    p = p + (step/6.0)*(kp1 + 2.0*kp2 + 2.0*kp3 + kp4);
  }
  kick = p*(pow(brho, 2.0));
}

void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double step, Mask& mask, Vector3D<> r, Vector3D<> p, std::vector<std::vector<double> >& trajectory){

  double alpha = 1.0/brho/beta;
  bool inside = true;
  Vector3D<> b; Vector3D<> b1; Vector3D<> b2; Vector3D<> b3;
  Vector3D<> kr1; Vector3D<> kp1; Vector3D<> r1; Vector3D<> p1;
  Vector3D<> kr2; Vector3D<> kp2; Vector3D<> r2; Vector3D<> p2;
  Vector3D<> kr3; Vector3D<> kp3; Vector3D<> r3; Vector3D<> p3;
  Vector3D<> kr4; Vector3D<> kp4;

  if (!trajectory.empty()) { trajectory.clear(); }
  std::vector<double> particle_pos;
  particle_pos.push_back(r.x); particle_pos.push_back(r.y); particle_pos.push_back(r.z);
  particle_pos.push_back(p.x); particle_pos.push_back(p.y); particle_pos.push_back(p.z);
  trajectory.push_back(particle_pos);
  particle_pos.clear();

  while (r.z < insertiondevice.z_max){

    inside = mask.is_inside(r); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b = insertiondevice.field(r);
    newton_lorentz_equation(alpha, r, p, b, kr1, kp1);
    r1 = r + (step/2.0)* kr1;
    p1 = p + (step/2.0)* kp1;

    inside = mask.is_inside(r1); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b1 = insertiondevice.field(r1);
    newton_lorentz_equation(alpha, r1, p1, b1, kr2, kp2);
    r2 = r + (step/2.0)* kr2;
    p2 = p + (step/2.0)* kp2;

    inside = mask.is_inside(r2); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b2 = insertiondevice.field(r2);
    newton_lorentz_equation(alpha, r2, p2, b2, kr3, kp3);
    r3 = r + step* kr3;
    p3 = p + step* kp3;

    inside = mask.is_inside(r3); if(!inside) { p.x = p.y = p.z = NAN; break; }
    b3 = insertiondevice.field(r3);
    newton_lorentz_equation(alpha, r3, p3, b3, kr4, kp4);

    r = r + (step/6.0)*(kr1 + 2.0*kr2 + 2.0*kr3 + kr4);
    p = p + (step/6.0)*(kp1 + 2.0*kp2 + 2.0*kp3 + kp4);

    particle_pos.push_back(r.x); particle_pos.push_back(r.y); particle_pos.push_back(r.z);
    particle_pos.push_back(p.x); particle_pos.push_back(p.y); particle_pos.push_back(p.z);
    trajectory.push_back(particle_pos);
    particle_pos.clear();

  }
}

void write_fieldmap_file(InsertionDevice& insertiondevice, std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  std::ofstream output_file(filename.c_str());

  if(!output_file){
    std::cout << "Can't open output file: " << filename <<  std::endl;
  } else {
    output_file << "Nome_do_Mapa:         " << "--" << std::endl;
    output_file << "Data_Hora:            " << "--" << std::endl;
    output_file << "Nome_do_Arquivo:      " << "--" << std::endl;
    output_file << "Numero_de_Imas:       " << "--" << std::endl;
    output_file << std::endl;
    output_file << "Nome_do_Ima:          " << "--" << std::endl;
    output_file << "Gap[mm]:              " << "--" << std::endl;
    output_file << "Gap_Controle[mm]:     " << "--" << std::endl;
    output_file << "Comprimento[mm]:      " << "--" << std::endl;
    output_file << "Corrente[A]:          " << "--" << std::endl;
    output_file << "Centro_Posicao_z[mm]: " << "--" << std::endl;
    output_file << "Centro_Posicao_x[mm]: " << "--" << std::endl;
    output_file << "Rotacao[graus]:       " << "--" << std::endl;
    output_file << std::endl;
    output_file << "X [mm]  Y[mm]  Z[mm] B x, B y, B z  [T]" << std::endl;
    output_file << "---------------------------------------------------------------------------------------------------------------" << std::endl;

    Vector3D<double> pos; Vector3D<double> b;
    for(int i = 0; i < z_vector.size(); i+=1){
      for (int j=0; j < y_vector.size(); j+=1){
        for(int k = 0; k < x_vector.size(); k+=1){
          pos.x = x_vector[k]; pos.y = y_vector[j]; pos.z = z_vector[i];
          b = insertiondevice.field(pos);
          output_file << std::scientific << std::showpos << pos.x*1000 << " ";
          output_file << std::scientific << std::showpos << pos.y*1000 << " ";
          output_file << std::scientific << std::showpos << pos.z*1000 << " ";
          output_file << std::scientific << std::showpos << b.x << " ";
          output_file << std::scientific << std::showpos << b.y << " ";
          output_file << std::scientific << std::showpos << b.z << std::endl;
        }
      }
    }

    std::cout << "Fieldmap saved in file: " << filename << std::endl << std::endl;
  }
}

void write_fieldmap_file(InsertionDevice& insertiondevice, std::string filename, std::vector<double> x_vector, double y_value, std::vector<double> z_vector){
  std::vector<double> y_vector; y_vector.push_back(y_value);
  write_fieldmap_file(insertiondevice, filename, x_vector, y_vector, z_vector);
}

void write_fieldmap_files(InsertionDevice& insertiondevice, std::string label, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  std::vector<double> y_value;
  std::ostringstream y_str;
  std::string filename_y;

  for(int i=0; i < y_vector.size(); i+=1){
    std::ostringstream ss;
    ss << (fabs(y_vector[i])*1000);
    std::string y_str(ss.str());

    if (y_vector[i] < 0){ filename_y = label  + "_yn" + y_str + ".txt"; }
    else if (y_vector[i] == 0){ filename_y = label + "_y0" + y_str + ".txt"; }
    else { filename_y = label + "_yp" + y_str + ".txt"; }

    y_value.push_back(y_vector[i]);
    write_fieldmap_file(insertiondevice, filename_y, x_vector, y_value, z_vector);
    y_value.clear();
  }

}

void save_kickmap(std::string filename, KickMap& kickmap){

  time_t now = time(0); char* date = ctime(&now);

  std::ofstream output_file(filename.c_str());
  if(!output_file) {
    std::cout << "Can't open output file: " << filename << std::endl << std::endl;
    return;
  }

  int nx = kickmap.x.size(); int ny = kickmap.y.size();

  output_file << "# KICKMAP" << std::endl;
  output_file << "# Author: Luana N. P. Vilela @ LNLS, Date: " << date;
  output_file << "# ID Length [m]" << std::endl;
  output_file << kickmap.physical_length << std::endl;
  output_file << "# Number of Horizontal Points" << std::endl;
  output_file << nx << std::endl;
  output_file << "# Number of Vertical Points" << std::endl;
  output_file << ny << std::endl;
  output_file << "# Horizontal KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";

  for(int j =0; j < nx; j+=1){
      output_file << std::scientific << std::showpos << kickmap.x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < ny; i+=1){
    output_file << std::scientific << std::showpos << kickmap.y[i] << " ";
    for(int j =0; j < nx; j+=1){
      if ( isnan(kickmap.kick_x[i][j]) ) { output_file << std::setw(13) << std::setfill(' ') << kickmap.kick_x[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << kickmap.kick_x[i][j] << " "; }
    }
    output_file << std::endl;
  }

  output_file << "# Vertical KickTable in T2m2" << std::endl;
  output_file << "START" << std::endl;
  output_file << std::setw(14) << " ";
  for(int j =0; j < nx; j+=1){
      output_file << std::scientific << std::showpos << kickmap.x[j] << " ";
  }
  output_file << std::endl;
  for(int i = 0; i < ny; i+=1){
    output_file << std::scientific << std::showpos << kickmap.y[i] << " ";
    for(int j =0; j < nx; j+=1){
      if ( isnan(kickmap.kick_y[i][j]) ){ output_file << std::setw(13) << std::setfill(' ') << kickmap.kick_y[i][j] << " "; }
      else { output_file << std::scientific << std::showpos << kickmap.kick_y[i][j] << " "; }
    }
    output_file << std::endl;
  }

  std::cout << "Kickmap saved in file: " << std::endl << filename << std::endl << std::endl;

}


void save_trajectories(std::string filename, std::vector<std::vector<std::vector<double> > >& trajectories){

  std::ofstream file(filename.c_str());
  if(!file) {
    std::cout << "Can't open output file: " << filename << std::endl << std::endl;
    return;
  }
  std::cout << std::endl << "Saving particle trajectories... " << std::endl;

  std::vector<std::vector<double> > trajectory;
  Vector3D<double> r;
  Vector3D<double> p;
  for (int i=0; i<trajectories.size(); i+=1){
    trajectory = trajectories[i];
    file << "START" << std::endl;
    for (int j=0; j<trajectory.size(); j+=1){
      r.x = trajectory[j][0]; r.y = trajectory[j][1]; r.z = trajectory[j][2];
      p.x = trajectory[j][3]; p.y = trajectory[j][3]; p.z = trajectory[j][5];
      file << r << " " << p << std::endl;
    }
  }

  std::cout << "Particle trajectories saved in file: " << std::endl;
  std::cout << filename << std::endl;

}

void save_trajectory(std::string filename, std::vector<std::vector<double> >& trajectory){

  std::ofstream file(filename.c_str());
  if(!file) {
    std::cout << "Can't open output file: " << filename << std::endl << std::endl;
    return;
  }
  std::cout << std::endl << "Saving particle trajectory... " << std::endl;

  Vector3D<double> r;
  Vector3D<double> p;
  for (int j=0; j<trajectory.size(); j+=1){
    r.x = trajectory[j][0]; r.y = trajectory[j][1]; r.z = trajectory[j][2];
    p.x = trajectory[j][3]; p.y = trajectory[j][3]; p.z = trajectory[j][5];
    file << r << " " << p << std::endl;
  }

  std::cout << "Particle trajectory saved in file: " << std::endl;
  std::cout << filename << std::endl;
}
