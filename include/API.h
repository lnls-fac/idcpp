#ifndef FIELDMAP_API_H
#define FIELDMAP_API_H

#include <vector>
#include "vector3d.hpp"
#include "fieldmap.h"

struct FieldMapException {
	enum type {
		success                       = 0,
		inconsistent_dimensions       = 1,
		out_of_range_x_min            = 2,
		out_of_range_x_max            = 3,
		out_of_range_z_min            = 4,
		out_of_range_z_max            = 5,
		file_not_found                = 6,
		out_of_range_y_min            = 7,
		out_of_range_y_max            = 8
	};
};

typedef std::vector<double> state_type;

#endif
