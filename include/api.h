#ifndef API_H
#define API_H

#include <string>
#include <vector>
#include <algorithm>
#include <alglib/linalg.h>
#include <alglib/interpolation.h>
#include "vector3d.hpp"
#include "matrix3d.h"


struct InsertionDeviceException {
    enum type {
        success = 0,
        inconsistent_dimensions = 1,
        file_not_found = 2,
        invalid_magnetization = 3,
				invalid_shape = 4,
				field_symmetry_error = 5,
    };
};

class Grid {
  public:
  	Grid() {};
  	Grid(int nx, int ny, double xmax, double ymax);
    Grid(int nx, int ny, double xmin, double xmax, double ymin, double ymax);
  	Grid(const Grid &obj);

  	int nx;
    int ny;
  	std::vector<double> x;
  	std::vector<double> y;

  private:
    double xmin, dx, xmax;
    double ymin, dy, ymax;

};


class Mask {
  public:
  	Mask() {};
    Mask(std::string shape, double width, double height);
    Mask(std::string filename);
  	Mask(const Mask &obj);

  	void load(std::string shape, double width, double height);
  	void load(std::string filename);
    bool is_inside(Vector3D<> position) const;


  private:
    double width;
    double height;
    std::vector<double> x;
    std::vector<double> y;
    std::string shape;
    alglib::spline1dinterpolant interpolant;

    void calc_interpolant();

};


class KickMap {
  public:
  	KickMap() {};
  	KickMap(double physical_length, std::vector<double> x, std::vector<double> y, std::vector<std::vector<double> > kick_x, std::vector<std::vector<double> > kicky);
    KickMap(const KickMap &obj);
    
  	double physical_length;
  	int nx;
  	int ny;
  	std::vector<double> x;
  	std::vector<double> y;
  	std::vector<std::vector<double> > kick_x;
  	std::vector<std::vector<double> > kick_y;
};


class Magnet {
  public:
    Magnet() {};
    virtual Vector3D<double> field(const Vector3D<double>& pos) const = 0;
    virtual double get_xmin() const = 0;
    virtual double get_xmax() const = 0;
    virtual double get_ymin() const = 0;
    virtual double get_ymax() const = 0;
    virtual double get_zmin() const = 0;
    virtual double get_zmax() const = 0;
    virtual double get_physical_length() const =0;
    std::vector<Vector3D<double> > field_vector(const std::vector<Vector3D<double> >& pos) const;
};


class FieldMap: public Magnet {
  public:

  	FieldMap(const std::string& fname_, size_t id_ = 0);
  	FieldMap(const FieldMap &obj);

  	size_t id;
  	size_t nx;
  	size_t nz;
    std::string fname;
    std::vector<double> x_grid;
  	std::vector<double> z_grid;

  	Vector3D<double> field(const Vector3D<double>& pos) const;
    size_t getid() const { return this->id; }
  	void   delete_data();
  	void   use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd);
    double get_xmin() const { return xmin; }
    double get_xmax() const { return xmax; }
    double get_ymin() const { return y; }
    double get_ymax() const { return y; }
    double get_zmin() const { return zmin; }
    double get_zmax() const { return zmax; }
    double get_physical_length() const { return physical_length; }

  private:
    double xmin, dx, xmax;
  	double zmin, dz, zmax;
  	double y;
  	double physical_length;
  	double *data;
  	alglib::spline2dinterpolant interpolant;

  	void calc_interpolant();
  	void read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions);

};


class FieldMap3D: public Magnet {
  public:

  	FieldMap3D() {};
  	FieldMap3D(const std::string& fname_, size_t id_ = 0);
  	FieldMap3D(const FieldMap3D &obj);

  	size_t id;
  	size_t nx;
  	size_t ny;
  	size_t nz;
  	std::string fname;
  	std::vector<double> x_grid;
  	std::vector<double> y_grid;
  	std::vector<double> z_grid;

  	Vector3D<double> field(const Vector3D<double>& pos) const;
  	size_t getid() const { return this->id; }
  	void   delete_data();
    double get_xmin() const { return xmin; }
    double get_xmax() const { return xmax; }
    double get_ymin() const { return ymin; }
    double get_ymax() const { return ymax; }
    double get_zmin() const { return zmin; }
    double get_zmax() const { return zmax; }
    double get_physical_length() const { return physical_length; }

  private:
    double xmin, dx, xmax;
  	double ymin, dy, ymax;
  	double zmin, dz, zmax;
  	double physical_length;
  	double *data;
  	std::vector<alglib::spline2dinterpolant> interpolants;

    void get_y_data(std::vector<double>& y_data, int y_index);
  	void calc_interpolants();
  	void read_fieldmap_from_file(const std::string& fname_, bool header, bool& consistent_dimensions);
  	Vector3D<double> field2D(const Vector3D<double>& pos, int interpolant_index) const;

};


class FieldMapContainer: public Magnet {
  public:

  	FieldMapContainer() {};
  	FieldMapContainer(FieldMap fieldmap);
  	FieldMapContainer(FieldMap3D fieldmap3D);
  	FieldMapContainer(std::vector<FieldMap> fieldmaps);
  	FieldMapContainer(std::string fieldmap_filename, bool fieldmap3D = false);
    FieldMapContainer(std::vector<std::string> fieldmap_filenames, bool fieldmap3D = false);

  	int nr_fieldmaps;

  	Vector3D<double> field(const Vector3D<double>& pos) const;
  	void use_field_symmetry(bool bx_odd, bool by_odd, bool bz_odd);
    double get_xmin() const { return xmin; }
    double get_xmax() const { return xmax; }
    double get_ymin() const { return ymin; }
    double get_ymax() const { return ymax; }
    double get_zmin() const { return zmin; }
    double get_zmax() const { return zmax; }
    double get_physical_length() const { return physical_length; }

  private:
    double xmin, xmax;
  	double ymin, ymax;
  	double zmin, zmax;
  	double physical_length;
    std::vector<FieldMap> fieldmaps;
  	FieldMap3D fieldmap3D;

    void set_attributes();
    void set_attributes3D();

};


class SubVolume {
  public:
    SubVolume( Vector3D<double> dim_ = 0.0, Vector3D<double> pos_ = 0.0, double str_ = 1.0) : pos(pos_), dim(dim_), str(str_) {};
    SubVolume(const SubVolume &obj);

  	Matrix3D<double> get_gmatrix(const Vector3D<double>& r) const;
    double           str;
    Vector3D<double> pos;
    Vector3D<double> dim;
};


class Block: public Magnet {
  public:
  	Block( Vector3D<double> mag_, Vector3D<double> dim_, Vector3D<double> pos_ = 0.0, double str_ = 1):
          mag(mag_) { subvolumes.push_back(SubVolume(dim_,pos_,str_)); }
    Block(const Block &obj);

  	Block&            add_subvolume( SubVolume sb) { subvolumes.push_back(sb); }
    Block&            add_subvolume( Vector3D<double> dim_, Vector3D<double> pos_, double str_ = -1) { subvolumes.push_back(SubVolume(dim_,pos_,str_)); }
    Vector3D<double>  field(const Vector3D<double>& r) const;
    Matrix3D<double>  get_gmatrix(const Vector3D<double>& r) const;
    Vector3D<double>  get_mag() const { return mag; }
    Vector3D<double>  get_pos() const { return subvolumes[0].pos; }
    Vector3D<double>  get_dim() const { return subvolumes[0].dim; }
    SubVolume         get_subvolume(unsigned int i) { return subvolumes[i]; }
    void 	 set_mag(Vector3D<double> mag_);
    void   set_pos(Vector3D<double> pos_);
    void   set_dim(Vector3D<double> dim_);
    void   set_subvolume(unsigned int i, SubVolume subvolume);
    virtual double get_xmin() const;
    double get_xmax() const;
    double get_ymin() const;
    double get_ymax() const;
    double get_zmin() const;
    double get_zmax() const;
    double get_physical_length() const;

  private:
    Vector3D<double>       mag;
    std::vector<SubVolume> subvolumes;

};


class BlockContainer : public Magnet {
  public:
    BlockContainer() {};
    Vector3D<double> field(const Vector3D<double>& r) const;
    BlockContainer& add_element(Block& b) { blocks.push_back(b); return *this; }
    Block& get_item(int i) { return blocks[i]; }
    BlockContainer& shift_pos( Vector3D<double> dr);
    int size() { return blocks.size(); }
    double get_xmin() const;
    double get_xmax() const;
    double get_ymin() const;
    double get_ymax() const;
    double get_zmin() const;
    double get_zmax() const;
    double get_physical_length() const;

  protected:
    std::vector<Block> blocks;

};


class HalbachCassette : public BlockContainer {
  public:
    HalbachCassette() {};
    HalbachCassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing=0.0,  int N=4);
  	HalbachCassette(const HalbachCassette &obj);

  	void gen_halbach_cassette(Block& genblock, const Matrix3D<double>& rot,  unsigned int nr_periods,  double spacing=0.0,  int N=4);
    HalbachCassette& 				set_xcenter(double x);
    HalbachCassette& 				set_zcenter(double z);
    HalbachCassette& 				set_ycenter(double y);
    HalbachCassette&        set_center_pos(Vector3D<double> pos);
    Vector3D<double>        get_center_pos() const;
    Vector3D<double> 				get_dim() const;
  	unsigned int						get_number_of_periods() {return nr_periods;}
  	double  								get_block_separation() {return spacing;}
  	int	 	   								get_number_of_blocks_per_period() {return N;}

  private:
  	unsigned int nr_periods;
  	double spacing;
  	int N;
};


class CassetteContainer: public Magnet {
  public:
    CassetteContainer() {};
    CassetteContainer(HalbachCassette& cassette);
    CassetteContainer(std::vector<HalbachCassette> cassettes);

    Vector3D<double> field(const Vector3D<double>& r) const;
    HalbachCassette& get_item(int i) { return cassettes[i]; }
    CassetteContainer& add_element(HalbachCassette& h) { cassettes.push_back(h); return *this; }
    int size() { return cassettes.size(); }
    double get_xmin() const;
    double get_xmax() const;
    double get_ymin() const;
    double get_ymax() const;
    double get_zmin() const;
    double get_zmax() const;
    double get_physical_length() const;

  private:
    std::vector<HalbachCassette> cassettes;
};


class EPU: public CassetteContainer {
  public:
    EPU() {};
    EPU(Block& genblock, unsigned int nr_periods, double magnetic_gap, double cassette_separation, double block_separation, double phase_csd, double phase_cie);
    void create_epu(Block& genblock, unsigned int nr_periods, double magnetic_gap, double cassette_separation, double block_separation, double phase_csd, double phase_cie);
    HalbachCassette& get_csd() {return get_item(0);}
    HalbachCassette& get_cse() {return get_item(1);}
    HalbachCassette& get_cid() {return get_item(2);}
    HalbachCassette& get_cie() {return get_item(3);}

  private:
    using CassetteContainer::add_element;
};

class DELTA: public CassetteContainer {
  public:
    DELTA() {};
    DELTA(Block& genblock, unsigned int nr_periods, double vertical_gap, double horizontal_gap, double block_separation, double phase_cd, double phase_ce);
    void create_delta(Block& genblock, unsigned int nr_periods, double vertical_gap, double horizontal_gap, double block_separation, double phase_cd, double phase_ce);
    HalbachCassette& get_cs() {return get_item(0);}
    HalbachCassette& get_cd() {return get_item(1);}
    HalbachCassette& get_ci() {return get_item(2);}
    HalbachCassette& get_ce() {return get_item(3);}

  private:
    using CassetteContainer::add_element;
};

void calc_brho(double energy, double& brho, double& beta);
void newton_lorentz_equation(double alpha, Vector3D<> r, Vector3D<> p,  Vector3D<> b, Vector3D<>& dr_ds, Vector3D<>& dp_ds);
void runge_kutta(Magnet& magnet, double brho, double beta, double step, Mask& mask, Vector3D<> r, Vector3D<> p, Vector3D<>& kick);
void runge_kutta(Magnet& magnet, double brho, double beta, double step, Mask& mask, Vector3D<> r, Vector3D<> p, std::vector<std::vector<double> >& trajectory);
void calc_kickmap(Magnet& magnet, Grid grid, Mask mask, double energy, double runge_kutta_step, KickMap& kickmap);
void calc_kickmap(Magnet& magnet, Grid grid, Mask mask, double energy, double runge_kutta_step, KickMap& kickmap, std::vector<std::vector<std::vector<double> > >& trajectories);

void write_fieldmap_file(Magnet& magnet, std::string filename, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);
void write_fieldmap_file(Magnet& magnet, std::string filename, std::vector<double> x_vector, double y_value, std::vector<double> z_vector);
void write_fieldmap_files(Magnet& magnet, std::string label, std::vector<double> x_vector, std::vector<double> y_vector, std::vector<double> z_vector);

void save_kickmap(std::string filename, KickMap& kickmap);
void save_trajectories(std::string filename, std::vector<std::vector<std::vector<double> > >& trajectories);
void save_trajectory(std::string filename, std::vector<std::vector<double> >& trajectory);


#endif
