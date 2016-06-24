#include <iostream>
#include <vector>
#include <cmath>
#include <api.h>

void get_field(std::vector<FieldMap>& fieldmaps, Vector3D<> r, Vector3D<>& b){
  if (fieldmaps.size() == 1){

    b = fieldmaps[0].field(r);

  } else {

    Vector3D<> f;
    std::vector<double> y_vector;
    std::vector<double> bx_vector;
    std::vector<double> by_vector;
    std::vector<double> bz_vector;

    for(int i=0; i < fieldmaps.size(); i+=1){
      f = fieldmaps[i].field(r);
      y_vector.push_back(fieldmaps[i].y);
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

    b.x = alglib::spline1dcalc(interpolant_x, r.y);
    b.y = alglib::spline1dcalc(interpolant_y, r.y);
    b.z = alglib::spline1dcalc(interpolant_z, r.y);
  }

}

void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds){
  dr_ds.x = p.x;
  dr_ds.y = p.y;
  dr_ds.z = p.z;
  dp_ds.x = - alpha * (p.y * b.z - p.z * b.y);
  dp_ds.y = - alpha * (p.z * b.x - p.x * b.z);
  dp_ds.z = - alpha * (p.x * b.y - p.y * b.x);
}

void runge_kutta(std::vector<FieldMap>& fieldmaps, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks){

  double alpha = 1.0/brho/beta;
  bool inside = true;
  Vector3D<> b; Vector3D<> b1; Vector3D<> b2; Vector3D<> b3;
  Vector3D<> kr1; Vector3D<> kp1; Vector3D<> r1; Vector3D<> p1;
  Vector3D<> kr2; Vector3D<> kp2; Vector3D<> r2; Vector3D<> p2;
  Vector3D<> kr3; Vector3D<> kp3; Vector3D<> r3; Vector3D<> p3;
  Vector3D<> kr4; Vector3D<> kp4;

  while (r.z < zmax){

    inside = mask.is_inside(r); if(!inside) { p.x = p.y = p.z = NAN; break; }
    get_field(fieldmaps, r, b);
    newton_lorentz_equation(alpha, r, p, b, kr1, kp1);
    r1 = r + (step/2.0)* kr1;
    p1 = p + (step/2.0)* kp1;

    inside = mask.is_inside(r1); if(!inside) { p.x = p.y = p.z = NAN; break; }
    get_field(fieldmaps, r1, b1);
    newton_lorentz_equation(alpha, r1, p1, b1, kr2, kp2);
    r2 = r + (step/2.0)* kr2;
    p2 = p + (step/2.0)* kp2;

    inside = mask.is_inside(r2); if(!inside) { p.x = p.y = p.z = NAN; break; }
    get_field(fieldmaps, r2, b2);
    newton_lorentz_equation(alpha, r2, p2, b2, kr3, kp3);
    r3 = r + step* kr3;
    p3 = p + step* kp3;

    inside = mask.is_inside(r3); if(!inside) { p.x = p.y = p.z = NAN; break; }
    get_field(fieldmaps, r3, b3);
    newton_lorentz_equation(alpha, r3, p3, b3, kr4, kp4);

    r = r + (step/6.0)*(kr1 + 2.0*kr2 + 2.0*kr3 + kr4);
    p = p + (step/6.0)*(kp1 + 2.0*kp2 + 2.0*kp3 + kp4);

  }

  kicks = p*(pow(brho, 2.0));

}
