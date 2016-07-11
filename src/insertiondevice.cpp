#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <api.h>

InsertionDevice::InsertionDevice(FieldMapContainer& fieldmap_container){
  this->fieldmaps = fieldmap_container;
  this->x_min = this->fieldmaps.x_min;
  this->x_max = this->fieldmaps.x_max;
  this->y_min = this->fieldmaps.y_min;
  this->y_max = this->fieldmaps.y_max;
  this->z_min = float(this->fieldmaps.z_min);
  this->z_max = float(this->fieldmaps.z_max);
  this->physical_length = this->fieldmaps.physical_length;
  this->type = 1;
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
  this->type = 1;
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
  this->type = 1;
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
  this->type = 2;
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
  this->type = 3;
}

Vector3D<double> InsertionDevice::field( Vector3D<double>& pos) {
  Vector3D<double> field;
  if (this->type == 1){
    field = this->fieldmaps.field(pos);
  } else if (this->type == 2){
    field = this->epu.field(pos);
  } else if (this->type == 3){
    field = this->delta.field(pos);
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

void InsertionDevice::write_fieldmap_file(std::string filename, std::vector<double> x_vector, double y, std::vector<double> z_vector){
  std::ofstream output_file(filename.c_str());

  if(! output_file){
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

    Vector3D<> pos; Vector3D<> field;
    for(int i = 0; i < z_vector.size(); i+=1){
      for(int k = 0; k < x_vector.size(); k+=1){
        pos.x = x_vector[k]; pos.y = y; pos.z = z_vector[i];
        field = this->field(pos);
        output_file << std::scientific << std::showpos << pos.x*1000 << " ";
        output_file << std::scientific << std::showpos << pos.y*1000 << " ";
        output_file << std::scientific << std::showpos << pos.z*1000 << " ";
        output_file << std::scientific << std::showpos << field.x << " ";
        output_file << std::scientific << std::showpos << field.y << " ";
        output_file << std::scientific << std::showpos << field.z << std::endl;
      }
    }

  std::cout << "Fieldmap saved in file: " << filename << std::endl;
  std::cout << std::endl;
  }
}

void InsertionDevice::write_fieldmap_files(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  double y;
  std::ostringstream y_str;
  std::string filename_y;

  for(int i=0; i < y_vector.size(); i+=1){
    y = y_vector[i];

    std::ostringstream ss; ss << (fabs(y)*1000);
    std::string y_str(ss.str());

    if (y < 0){ filename_y = filename  + "_yn" + y_str + ".txt"; }
    else if (y == 0){ filename_y = filename + "_y0" + y_str + ".txt"; }
    else { filename_y = filename + "_yp" + y_str + ".txt"; }

    this->write_fieldmap_file(filename_y, x_vector, y, z_vector);
  }
}

void InsertionDevice::write_fieldmap3D_file(std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector){
  std::ofstream output_file(filename.c_str());

  if(! output_file){
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

    Vector3D<> pos; Vector3D<> field;
    for(int i = 0; i < z_vector.size(); i+=1){
      for (int j=0; j < y_vector.size(); j+=1){
        for(int k = 0; k < x_vector.size(); k+=1){
          pos.x = x_vector[k]; pos.y = y_vector[j]; pos.z = z_vector[i];
          field = this->field(pos);
          output_file << std::scientific << std::showpos << pos.x*1000 << " ";
          output_file << std::scientific << std::showpos << pos.y*1000 << " ";
          output_file << std::scientific << std::showpos << pos.z*1000 << " ";
          output_file << std::scientific << std::showpos << field.x << " ";
          output_file << std::scientific << std::showpos << field.y << " ";
          output_file << std::scientific << std::showpos << field.z << std::endl;
        }
      }
    }

    std::cout << "Fieldmap saved in file: " << filename << std::endl;
    std::cout << std::endl;
  }
}

void InsertionDevice::calc_kickmap(Grid grid, Mask mask, double energy, double runge_kutta_step){

  double beta; double brho;
  calc_brho(energy, brho, beta);

  double z_min = this->z_min;
  double z_max = this->z_max;

  int count = 0;
  int size = grid.nx * grid.ny;
  std::vector<double> kick_x_vector;
  std::vector<double> kick_y_vector;
  Vector3D<> r(0.0, 0.0, z_min+1);
  Vector3D<> p(0.0, 0.0, 1.0);
  Vector3D<> kick;

  std::cout << std::endl;
  std::cout << "Calculating kickmap..." << std::endl;

  std::vector<std::vector<double> > kick_x;
	std::vector<std::vector<double> > kick_y;

  for(int i = 0; i < grid.ny; i+=1){
    kick_x_vector.clear();
    kick_y_vector.clear();

    for(int j =0; j < grid.nx; j+=1){
      r.x = grid.x[j];
      r.y = grid.y[i];
      runge_kutta(*this, brho, beta, z_max, runge_kutta_step, mask, r, p, kick);
      kick_x_vector.push_back(kick.x);
      kick_y_vector.push_back(kick.y);

      count += 1;
      if (count%10 == 0) {
        std::cout << std::setw(6) << std::setprecision(4) << std::setfill(' ') << 100.0*(double(count)/double(size)) << '%' << '\r' << std::flush;
      }

    }
    kick_x.push_back(kick_x_vector);
    kick_y.push_back(kick_y_vector);
  }

  KickMap kickmap(this->physical_length, grid.x, grid.y, kick_x, kick_y);
  this->kickmap = kickmap;

}

void InsertionDevice::write_kickmap_file(std::string filename){
  this->kickmap.write_to_file(filename);
}
