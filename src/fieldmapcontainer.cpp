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
    this->set_attributes3D();

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
        std::cout << "Load " << filename << ":  y = " << fieldmap.get_ymin() << std::endl;
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
      std::cout << "Load " << fieldmap_filename << ":  y = " << fieldmap.get_ymin() << std::endl;
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
    this->set_attributes3D();
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
  this->set_attributes3D();
}

void FieldMapContainer::set_attributes(){
  this->nr_fieldmaps = this->fieldmaps.size();
  std::vector<double> xmin_vector;
  std::vector<double> xmax_vector;
  std::vector<double> y_vector;
  std::vector<double> zmin_vector;
  std::vector<double> zmax_vector;
  for(int i=0; i < this->nr_fieldmaps; i+=1) {
    xmin_vector.push_back(this->fieldmaps[i].get_xmin());
    xmax_vector.push_back(this->fieldmaps[i].get_xmax());
    y_vector.push_back(this->fieldmaps[i].get_ymin());
    zmin_vector.push_back(this->fieldmaps[i].get_zmin());
    zmax_vector.push_back(this->fieldmaps[i].get_zmax());
  }
  this->xmin = *std::min_element(xmin_vector.begin(), xmin_vector.end());
  this->xmax = *std::max_element(xmax_vector.begin(), xmax_vector.end());
  this->ymin = *std::min_element(y_vector.begin(), y_vector.end());
  this->ymax = *std::max_element(y_vector.begin(), y_vector.end());
  this->zmin = *std::min_element(zmin_vector.begin(), zmin_vector.end());
  this->zmax = *std::max_element(zmax_vector.begin(), zmax_vector.end());
  this->physical_length = this->fieldmaps[0].get_physical_length();
}

void FieldMapContainer::set_attributes3D(){
  this->nr_fieldmaps = 0;
  this->xmin = this->fieldmap3D.get_xmin();
  this->xmax = this->fieldmap3D.get_xmax();
  this->ymin = this->fieldmap3D.get_ymin();
  this->ymax = this->fieldmap3D.get_ymax();
  this->zmin = this->fieldmap3D.get_zmin();
  this->zmax = this->fieldmap3D.get_zmax();
  this->physical_length = this->fieldmap3D.get_physical_length();
}

void FieldMapContainer::use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd){
  std::vector<FieldMap> new_fieldmaps;

  if (this->nr_fieldmaps == 0){ throw InsertionDeviceException::field_symmetry_error; }

  int negative_count = 0; int positive_count = 0; int zero_count = 0;
  for (int i=0; i < this->nr_fieldmaps; i+=1){
    if (this->fieldmaps[i].get_ymin() < 0) { negative_count +=1;}
    else if (this->fieldmaps[i].get_ymin() > 0) { positive_count +=1;}
    else if (this->fieldmaps[i].get_ymin() == 0) { zero_count +=1;}
  }
  if (((negative_count+zero_count)!= this->nr_fieldmaps) && ((positive_count+zero_count)!= this->nr_fieldmaps)){
    throw InsertionDeviceException::field_symmetry_error;
  }

  for (int i = (this->nr_fieldmaps -1); i >= 0; i-=1){
    if (this->fieldmaps[i].get_ymin() != 0){
      FieldMap neg_fieldmap(this->fieldmaps[i]);
      neg_fieldmap.use_field_symmetry(bx_odd, by_odd, bz_odd);
      new_fieldmaps.push_back(neg_fieldmap);
      std::cout << "Load " << neg_fieldmap.fname << ":  y = " << neg_fieldmap.get_ymin() << std::endl;
    }
  }

  for (int i=0; i < this->nr_fieldmaps; i+=1){
    new_fieldmaps.push_back(this->fieldmaps[i]);
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
      y_vector.push_back(this->fieldmaps[i].get_ymin());
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
