#include <string>
#include <iostream>

#include "fieldmap.h"


std::string fname = "/home/fac_files/misc/c++/insertion_devices/fmap.txt";
Vector3D<> calc_field(Vector3D<> r);

int main() {

  FieldMap fmap(fname);

  Vector3D<> r(0,0,0);
  Vector3D<> b = fmap.field(r);

  std::cout << b.y << std::endl;

  b = calc_field(r);

}


Vector3D<> calc_field(Vector3D<> r) {
  Vector3D<> b;
  // calc
  return b;
}
