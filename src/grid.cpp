#include <vector>
#include <api.h>

Grid::Grid(int nx, int ny, double xmax, double ymax){

  this->nx = nx;
  this->ny = ny;
  this->xmax = xmax;
  this->ymax = ymax;
  this->xmin = (this->nx > 1) ? - xmax : xmax;
  this->ymin = (this->ny > 1) ? - ymax : ymax;

  this->dx = (this->nx > 1) ? (this->xmax - this->xmin)/(this->nx - 1) : 0.0;
  this->dy = (this->ny > 1) ? (this->ymax - this->ymin)/(this->ny - 1) : 0.0;

  this->x.clear();
  this->y.clear();

  if (this->nx > 1){
    for(int i = 0; i < this->nx; i+=1){
      this->x.push_back( this->xmin + double(i*this->dx) );
    }
  } else this->x.push_back(this->xmax);
  if (this->ny > 1){
    for(int j = 0; j < this->ny; j+=1){
      this->y.push_back( this->ymin + double(j*this->dy) );
    }
  } else this->y.push_back(this->ymax);

}


Grid::Grid(int nx, int ny, double xmin, double xmax, double ymin, double ymax){

  if (nx == 1 && xmin != xmax){throw InsertionDeviceException::inconsistent_dimensions; }
  if (ny == 1 && ymin != ymax){throw InsertionDeviceException::inconsistent_dimensions; }

  this->nx = nx;
  this->ny = ny;
  this->xmax = xmax;
  this->ymax = ymax;
  this->xmin = xmin;
  this->ymin = ymin;

  this->dx = (this->nx > 1) ? (this->xmax - this->xmin)/(this->nx - 1) : 0.0;
  this->dy = (this->ny > 1) ? (this->ymax - this->ymin)/(this->ny - 1) : 0.0;

  this->x.clear();
  this->y.clear();

  if (this->nx > 1){
    for(int i = 0; i < this->nx; i+=1){
      this->x.push_back( this->xmin + double(i*this->dx) );
    }
  } else this->x.push_back(this->xmax);
  if (this->ny > 1){
    for(int j = 0; j < this->ny; j+=1){
      this->y.push_back( this->ymin + double(j*this->dy) );
    }
  } else this->y.push_back(this->ymax);

}

Grid::Grid(const Grid &obj){
  this->nx = obj.nx;        this->ny = obj.ny;
  this->xmax = obj.xmax;    this->ymax = obj.ymax;
  this->xmin = obj.xmin;    this->ymin = obj.ymin;
  this->dx = obj.dx;        this->dy = obj.dy;
  this->x = obj.x;          this->y = obj.y;
}
