#include <vector>
#include <api.h>

struct GridException { enum type { inconsistent_dimensions = 0,}; };

Grid::Grid(int nx, int ny, double x_max, double y_max){
  //symmetrical grid

  this->nx = nx;
  this->ny = ny;
  this->x_max = x_max;
  this->y_max = y_max;
  this->x_min = (this->nx > 1) ? - x_max : x_max;
  this->y_min = (this->ny > 1) ? - y_max : y_max;

  this->dx = (this->nx > 1) ? (this->x_max - this->x_min)/(this->nx - 1) : 0;
  this->dy = (this->ny > 1) ? (this->y_max - this->y_min)/(this->ny - 1) : 0;

  this->x.clear();
  this->y.clear();

  if (this->nx > 1){
    for(int i = 0; i < this->nx; i+=1){
      this->x.push_back( this->x_min + double(i*this->dx) );
    }
  } else this->x.push_back(this->x_max);
  if (this->ny > 1){
    for(int j = 0; j < this->ny; j+=1){
      this->y.push_back( this->y_min + double(j*this->dy) );
    }
  } else this->y.push_back(this->y_max);

}


Grid::Grid(int nx, int ny, double x_min, double x_max, double y_min, double y_max){

  if (nx == 1 && x_min != x_max){throw GridException::inconsistent_dimensions; }
  if (ny == 1 && y_min != y_max){throw GridException::inconsistent_dimensions; }

  this->nx = nx;
  this->ny = ny;
  this->x_max = x_max;
  this->y_max = y_max;
  this->x_min = x_min;
  this->y_min = y_min;

  this->dx = (this->nx > 1) ? (this->x_max - this->x_min)/(this->nx - 1) : 0;
  this->dy = (this->ny > 1) ? (this->y_max - this->y_min)/(this->ny - 1) : 0;

  this->x.clear();
  this->y.clear();

  if (this->nx > 1){
    for(int i = 0; i < this->nx; i+=1){
      this->x.push_back( this->x_min + double(i*this->dx) );
    }
  } else this->x.push_back(this->x_max);
  if (this->ny > 1){
    for(int j = 0; j < this->ny; j+=1){
      this->y.push_back( this->y_min + double(j*this->dy) );
    }
  } else this->y.push_back(this->y_max);

}

Grid::Grid(const Grid &obj){
  this->nx = obj.nx;        this->ny = obj.ny;
  this->x_max = obj.x_max;  this->y_max = obj.y_max;
  this->x_min = obj.x_min;  this->y_min = obj.y_min;
  this->dx = obj.dx;        this->dy = obj.dy;
  this->x = obj.x;          this->y = obj.y;
}
