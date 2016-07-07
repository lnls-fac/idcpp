%module insertion_devices

%{
    #include "api.h"
    #include "vector3d.h"
%}

%include "std_string.i"
%include "std_vector.i"


%include "api.h"
%include "vector3d.h"

%template(CppVector3D) Vector3D<double>;

namespace std {
      %template(CppStringVector) vector<string>;
      %template(CppDoubleVector) vector<double>;
      %template(CppVectorVector3D) vector<Vector3D<double> >;
}
