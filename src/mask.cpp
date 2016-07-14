#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <api.h>

enum Shape{ SQUARE, RECTANGLE, DIAMOND, CIRCLE, ELLIPSE, TABLE, NONE};
std::map<std::string, Shape> shapes;

void defined_shapes(){
  shapes["RECTANGLE"] = RECTANGLE;
  shapes["DIAMOND"]   = DIAMOND;
  shapes["ELLIPSE"]   = ELLIPSE;
  shapes["TABLE"]     = TABLE;
  shapes["NONE"]      = NONE;
}

Mask::Mask(std::string shape, double width=0.0, double height=0.0){ this->load(shape, width, height); }

Mask::Mask(std::string filename){ this->load(filename); }

Mask::Mask(const Mask &obj){
  this->shape = obj.shape;
  this->width = obj.width;
  this->height = obj.height;
  this->x = obj.x;
  this->y = obj.y;
  this->interpolant = obj.interpolant;
}

void Mask::load(std::string shape, double width, double height){
  std::transform(shape.begin(), shape.end(), shape.begin(), ::toupper);
  defined_shapes();

  if (shapes.find(shape) == shapes.end()) throw InsertionDeviceException::invalid_shape;

  this->shape = shape;
  this->width = width;
  this->height = height;
}

void Mask::load(std::string filename){
  std::ifstream file(filename.c_str());
  if (!file.is_open()) throw InsertionDeviceException::file_not_found;

  std::vector<double> x; std::vector<double> y;
  double x_value; double y_value;
  std::string line; std::getline(file, line); std::getline(file, line);

  while (true) {
    file >> x_value >> y_value;
    if (file.eof()) break;
    if (x.size() >= 1){
      // Avoid alglib error
      if (x.back() == (x_value/1000.0)){ x_value += 1e-10; }
      if (y.back() == (y_value/1000.0)){ y_value += 1e-10; }
    }
    x.push_back(double(x_value)/1000.0);
    y.push_back(double(y_value)/1000.0);
  }
  file.close();

  this->shape = "TABLE";
  this->x = x;
  this->y = y;
  this->calc_interpolant();
  defined_shapes();
}

void Mask::calc_interpolant(){
  alglib::real_1d_array x_array;
  alglib::real_1d_array y_array;

  x_array.setcontent(this->x.size(), &this->x[0]);
  y_array.setcontent(this->y.size(), &this->y[0]);

  alglib::spline1dinterpolant interpolant;
  alglib::spline1dbuildlinear(x_array, y_array, interpolant);
  this->interpolant = interpolant;
}

bool Mask::is_inside(Vector3D<> position) const{
  bool inside;
  double pos_x = fabs(position.x);
  double pos_y = fabs(position.y);

  switch(shapes[this->shape]){
    case RECTANGLE:
      inside = ((pos_x <= (this->width/2.0)) && (pos_y <= (this->height/2.0))) ? true : false;
      break;
    case DIAMOND:
      inside = (pos_y <= (this->height/2.0) - (this->height/this->width)*pos_x) ? true : false;
      break;
    case ELLIPSE:
      inside = (pow((pos_x/(this->width/2.0)), 2.0) + pow((pos_y/(this->height/2.0)), 2.0) <= 1.0) ? true : false;
      break;
    case NONE:
      inside = true;
      break;
    case TABLE:
      double x_min = *std::min_element(this->x.begin(), this->x.end());
      double x_max = *std::max_element(this->x.begin(), this->x.end());
      inside = ((pos_x >= x_min) && (pos_x <= x_max) && \
               (pos_y <= alglib::spline1dcalc(this->interpolant, pos_x))) ? true : false;
      break;
  }

  return inside;
}
