%module showcase
/*SWIG has a library for doing the string mappings.*/
%include "std_string.i"
%{

/* Includes the header in the wrapper code */
#include "libshowcase.h"
%}
 
/* Parse the header file to generate wrappers */
%include "libshowcase.h"