#include <api.h>

KickMap::KickMap(double physical_length, std::vector<double> x, std::vector<double> y, std::vector<std::vector<double> > kick_x, std::vector<std::vector<double> > kick_y){
  this->physical_length = physical_length;
  this->nx = x.size();
  this->ny = y.size();
  this->x = x;
  this->y = y;
  this->kick_x = kick_x;
  this->kick_y = kick_y;
}


KickMap::KickMap(const KickMap &obj){
  this->physical_length = obj.physical_length;
  this->nx = obj.nx;
  this->ny = obj.ny;
  this->x = obj.x;
  this->y = obj.y;
  this->kick_x = obj.kick_x;
  this->kick_y = obj.kick_y;
}
