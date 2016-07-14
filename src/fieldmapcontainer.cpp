#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <api.h>

FieldMapContainer::FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool fieldmap3D){
  std::string filename;

  if ((fieldmap3D) && (fieldmap_filenames.size() != 1)){

    throw InsertionDeviceException::inconsistent_dimensions;

  } else if ((fieldmap3D) && (fieldmap_filenames.size() == 1)){

    std::cout << "Loading fieldmap file..." << std::endl;
    try{
      FieldMap3D fieldmap3D(filename.c_str());
      this->fieldmap3D = fieldmap3D;
      std::cout << "Load " << filename << std::endl;
    } catch(InsertionDeviceException::type e){
      if (e == 2) { std::cout << "File not found: " << filename << std::endl << std::endl; }
      throw InsertionDeviceException::file_not_found;
    }
    this->nr_fieldmaps = 0;
    this->x_min = this->fieldmap3D.x_min;
    this->x_max = this->fieldmap3D.x_max;
    this->y_min = this->fieldmap3D.y_min;
    this->y_max = this->fieldmap3D.y_max;
    this->z_min = this->fieldmap3D.z_min;
    this->z_max = this->fieldmap3D.z_max;
    this->physical_length = this->fieldmap3D.physical_length;

  } else {

    int nr_fieldmaps = 0;
    int nr_files = fieldmap_filenames.size();

    std::cout << "Loading fieldmap files..." << std::endl;
    for(int i=0; i < nr_files; i+=1){
      try {
        filename = fieldmap_filenames.back();
        fieldmap_filenames.pop_back();
        FieldMap fieldmap(filename.c_str());
        this->fieldmaps.push_back(fieldmap);
        std::cout << "Load " << filename << ":  y = " << fieldmap.y << std::endl;
        nr_fieldmaps +=1;
      } catch(InsertionDeviceException::type e){
        if (e == 1) { std::cout << "File not found: " << filename << std::endl << std::endl; }
        throw InsertionDeviceException::file_not_found;
      }
    }
    this->set_attributes();
  }
}

FieldMapContainer::FieldMapContainer(std::string fieldmap_filename, bool fieldmap3D){
  std::cout << "Loading fieldmap file..." << std::endl;

  if (!fieldmap3D){
    try{
      FieldMap fieldmap(fieldmap_filename.c_str());
      this->fieldmaps.push_back(fieldmap);
      std::cout << "Load " << fieldmap_filename << ":  y = " << fieldmap.y << std::endl;
    } catch(InsertionDeviceException::type e){
      if (e == 2) { std::cout << "File not found: " << fieldmap_filename << std::endl << std::endl; }
      throw InsertionDeviceException::file_not_found;
    }
    this->set_attributes();
  } else {
    try{
      FieldMap3D fieldmap3D(fieldmap_filename.c_str());
      this->fieldmap3D = fieldmap3D;
      std::cout << "Load " << fieldmap_filename << std::endl;
    } catch(InsertionDeviceException::type e){
      if (e == 2) { std::cout << "File not found: " << fieldmap_filename << std::endl << std::endl; }
      throw InsertionDeviceException::file_not_found;
    }
    this->nr_fieldmaps = 0;
    this->x_min = this->fieldmap3D.x_min;
    this->x_max = this->fieldmap3D.x_max;
    this->y_min = this->fieldmap3D.y_min;
    this->y_max = this->fieldmap3D.y_max;
    this->z_min = this->fieldmap3D.z_min;
    this->z_max = this->fieldmap3D.z_max;
    this->physical_length = this->fieldmap3D.physical_length;
  }

}

FieldMapContainer::FieldMapContainer(std::vector<FieldMap> fieldmaps){
  this->fieldmaps = fieldmaps;
  this->set_attributes();
}

FieldMapContainer::FieldMapContainer(FieldMap fieldmap){
  this->fieldmaps.push_back(fieldmap);
  this->set_attributes();
}

FieldMapContainer::FieldMapContainer(FieldMap3D fieldmap3D){
  this->fieldmap3D = fieldmap3D;
  this->nr_fieldmaps = 0;
  this->x_min = this->fieldmap3D.x_min;
  this->x_max = this->fieldmap3D.x_max;
  this->y_min = this->fieldmap3D.y_min;
  this->y_max = this->fieldmap3D.y_max;
  this->z_min = this->fieldmap3D.z_min;
  this->z_max = this->fieldmap3D.z_max;
  this->physical_length = this->fieldmap3D.physical_length;
}

void FieldMapContainer::set_attributes(){
  this->nr_fieldmaps = this->fieldmaps.size();
  std::vector<double> x_min_vector;
  std::vector<double> x_max_vector;
  std::vector<double> y_vector;
  std::vector<double> z_min_vector;
  std::vector<double> z_max_vector;
  for(int i=0; i < this->nr_fieldmaps; i+=1) {
    x_min_vector.push_back(this->fieldmaps[i].x_min);
    x_max_vector.push_back(this->fieldmaps[i].x_max);
    y_vector.push_back(this->fieldmaps[i].y);
    z_min_vector.push_back(this->fieldmaps[i].z_min);
    z_max_vector.push_back(this->fieldmaps[i].z_max);
  }
  this->x_min = *std::min_element(x_min_vector.begin(), x_min_vector.end());
  this->x_max = *std::max_element(x_max_vector.begin(), x_max_vector.end());
  this->y_min = *std::min_element(y_vector.begin(), y_vector.end());
  this->y_max = *std::max_element(y_vector.begin(), y_vector.end());
  this->z_min = *std::min_element(z_min_vector.begin(), z_min_vector.end());
  this->z_max = *std::max_element(z_max_vector.begin(), z_max_vector.end());
  this->physical_length = this->fieldmaps[0].physical_length;
}

void FieldMapContainer::use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd){
  std::vector<FieldMap> new_fieldmaps;

  if (this->nr_fieldmaps == 0){ throw InsertionDeviceException::field_symmetry_error; }

  int negative_count = 0; int positive_count = 0; int zero_count = 0;
  for (int i=0; i < this->nr_fieldmaps; i+=1){
    if (this->fieldmaps[i].y < 0) { negative_count +=1;}
    else if (this->fieldmaps[i].y > 0) { positive_count +=1;}
    else if (this->fieldmaps[i].y == 0) { zero_count +=1;}
  }
  if (((negative_count+zero_count)!= this->nr_fieldmaps) && ((positive_count+zero_count)!= this->nr_fieldmaps)){
    throw InsertionDeviceException::field_symmetry_error;
  }

  for (int i = (this->nr_fieldmaps -1); i >= 0; i-=1){
    if (this->fieldmaps[i].y != 0){
      FieldMap neg_fieldmap(this->fieldmaps[i]);
      neg_fieldmap.use_field_symmetry(bx_odd, by_odd, bz_odd);
      new_fieldmaps.push_back(neg_fieldmap);
      std::cout << "Load " << neg_fieldmap.fname << ":  y = " << neg_fieldmap.y << std::endl;
    }
  }

  for (int i=0; i < this->nr_fieldmaps; i+=1){
    new_fieldmaps.push_back(this->fieldmaps[i]);
  }

  //Teste!
  for (int k=0; k < new_fieldmaps.size(); k+=1){
    std::cout << new_fieldmaps[k].y << std::endl;
  }

  this->fieldmaps = new_fieldmaps;
  this->set_attributes();

}

Vector3D<double> FieldMapContainer::field(const Vector3D<>& pos) const{

  Vector3D<> field;

  if (this->nr_fieldmaps == 0){

    field = this->fieldmap3D.field(pos);

  } else if (this->nr_fieldmaps == 1){

    field = this->fieldmaps[0].field(pos);

  } else {

    Vector3D<> f;
    std::vector<double> y_vector;
    std::vector<double> bx_vector;
    std::vector<double> by_vector;
    std::vector<double> bz_vector;

    for(int i=0; i < this->fieldmaps.size(); i+=1){
      f = this->fieldmaps[i].field(pos);
      y_vector.push_back(this->fieldmaps[i].y);
      bx_vector.push_back(f.x);
      by_vector.push_back(f.y);
      bz_vector.push_back(f.z);
    }

    alglib::real_1d_array y_array;
    alglib::real_1d_array bx_array;
    alglib::real_1d_array by_array;
    alglib::real_1d_array bz_array;
    y_array.setcontent(y_vector.size(), &y_vector[0]);
    bx_array.setcontent(bx_vector.size(), &bx_vector[0]);
    by_array.setcontent(by_vector.size(), &by_vector[0]);
    bz_array.setcontent(bz_vector.size(), &bz_vector[0]);

    alglib::spline1dinterpolant interpolant_x;
    alglib::spline1dinterpolant interpolant_y;
    alglib::spline1dinterpolant interpolant_z;
    alglib::spline1dbuildcubic(y_array, bx_array, interpolant_x);
    alglib::spline1dbuildcubic(y_array, by_array, interpolant_y);
    alglib::spline1dbuildcubic(y_array, bz_array, interpolant_z);

    field.x = alglib::spline1dcalc(interpolant_x, pos.y);
    field.y = alglib::spline1dcalc(interpolant_y, pos.y);
    field.z = alglib::spline1dcalc(interpolant_z, pos.y);
  }

  return field;

}

std::vector<Vector3D<double> > FieldMapContainer::field(const std::vector<Vector3D<> >& pos) const{

  std::vector<Vector3D<> > field;
  for (int i=0; i < pos.size(); i+=1){
    field.push_back(this->field(pos[i]));
  }
  return field;

}
