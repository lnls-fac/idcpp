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
