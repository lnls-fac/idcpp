%module insertion_devices

%{
    #include "api.h"
%}

%include "std_string.i"
%include "std_vector.i"


%include "api.h"

namespace std {
      %template(CppStringVector) vector<string>;
      %template(CppDoubleVector) vector<double>;
//    %template(CppPVValuePairVector) vector<PVValuePair>;
}
