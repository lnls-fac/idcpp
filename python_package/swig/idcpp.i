%module idcpp

%{
    #include "api.h"
    #include "halbachcassette.h"
    #include "vector3d.hpp"
    #include "matrix3d.h"
%}

%include "std_string.i"
%include "std_vector.i"

%include "api.h"
%include "halbachcassette.h"
%include "vector3d.hpp"
%include "matrix3d.h"

%template(CppVector3D) Vector3D<double>;
%template(CppMatrix3D) Matrix3D<double>;

namespace std {
      %template(CppStringVector) vector<string>;
      %template(CppDoubleVector) vector<double>;
      %template(CppVectorVector3D) vector<Vector3D<double> >;
      %template(CppDoubleVectorVector) vector<vector<double> >;
}
