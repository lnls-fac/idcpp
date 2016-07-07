#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <api.h>

FieldMapContainer::FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool use_field_simmetry){

  std::string filename;
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

      if (use_field_simmetry && (fieldmap.y > 0)) {
        FieldMap neg_fieldmap(fieldmap);
        neg_fieldmap.change_y_sign();
        this->fieldmaps.push_back(neg_fieldmap);
        std::cout << "Load " << filename << ":  y = " << neg_fieldmap.y << std::endl;
        nr_fieldmaps +=1;
      }

    }
    catch(FieldMapException::type e){
      if (e == 1) { std::cout << "File not found: " << filename << std::endl << std::endl; }
      throw;
    }
  }
  this->set_attributes();
}

FieldMapContainer::FieldMapContainer(std::vector<FieldMap> fieldmaps){
  this->fieldmaps = fieldmaps;
  this->set_attributes();
}

FieldMapContainer::FieldMapContainer(FieldMap fieldmap){
  this->fieldmaps.push_back(fieldmap);
  this->set_attributes();
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

Vector3D<double> FieldMapContainer::field(const Vector3D<>& pos) const{

  Vector3D<> field;

  if (this->nr_fieldmaps == 1){

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
