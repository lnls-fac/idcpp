#include <iostream>
#include <vector>
#include <cmath>
#include <api.h>

void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds){

  dr_ds.x = p.x;
  dr_ds.y = p.y;
  dr_ds.z = p.z;
  dp_ds.x = -alpha * (p.y * b.z - p.z * b.y);
  dp_ds.y = -alpha * (p.z * b.x - p.x * b.z);
  dp_ds.z = -alpha * (p.x * b.y - p.y * b.x);

}

void runge_kutta(InsertionDevice& insertiondevice, double brho, double beta, double zmax, double step, Mask mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kicks){

  double alpha = 1.0/brho/beta;
  bool inside = true;
  Vector3D<> b; Vector3D<> b1; Vector3D<> b2; Vector3D<> b3;
  Vector3D<> kr1; Vector3D<> kp1; Vector3D<> r1; Vector3D<> p1;
  Vector3D<> kr2; Vector3D<> kp2; Vector3D<> r2; Vector3D<> p2;
  Vector3D<> kr3; Vector3D<> kp3; Vector3D<> r3; Vector3D<> p3;
  Vector3D<> kr4; Vector3D<> kp4;

  while (r.z < zmax){

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

  kicks = p*(pow(brho, 2.0));

}
