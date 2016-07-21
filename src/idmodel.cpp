#include <api.h>

EPU::EPU(Block& genblock, unsigned int nr_periods, double magnetic_gap, double cassette_separation, double block_separation, double phase_csd, double phase_cie){
  this->create_epu(genblock, nr_periods, magnetic_gap, cassette_separation, block_separation, phase_csd, phase_cie);
}

void EPU::create_epu(Block& genblock, unsigned int nr_periods, double magnetic_gap, double cassette_separation, double block_separation, double phase_csd, double phase_cie){

  HalbachCassette csd; HalbachCassette cse; HalbachCassette cid; HalbachCassette cie;
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

  this->add_element(csd);
  this->add_element(cse);
  this->add_element(cid);
  this->add_element(cie);
}


DELTA::DELTA(Block& genblock, unsigned int nr_periods, double vertical_gap, double horizontal_gap, double block_separation, double phase_cd, double phase_ce){
  this->create_delta(genblock, nr_periods, vertical_gap, horizontal_gap, block_separation, phase_cd, phase_ce);
}

void DELTA::create_delta(Block& genblock, unsigned int nr_periods, double vertical_gap, double horizontal_gap, double block_separation, double phase_cd, double phase_ce){

  HalbachCassette cs; HalbachCassette cd; HalbachCassette ci; HalbachCassette ce;
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

  this->add_element(cs);
  this->add_element(cd);
  this->add_element(ci);
  this->add_element(ce);
}
