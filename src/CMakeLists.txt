INCLUDE(CheckFunctionExists)

add_subdirectory(libmarqov)

# configuration for the simple ising test case
add_executable(mysimpleising mysimpleising.cpp libmarqov/marqov.cpp)
target_include_directories(mysimpleising PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(mysimpleising libmarqov)
set_target_properties(mysimpleising PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

# configuration for the Heisenberg test case with threading
add_executable(mysimpleheisenbergthreaded mysimpleheisenbergthreaded.cpp libmarqov/marqov.cpp)
target_include_directories(mysimpleheisenbergthreaded PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(mysimpleheisenbergthreaded libmarqov)
set_target_properties(mysimpleheisenbergthreaded PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

# configuration of the main parts.

SET(MARQOVTARGETS MARQOV mpiMARQOV)

foreach(target ${MARQOVTARGETS})
add_executable (${target} libmarqov/util/registry.cpp main.cpp)
target_link_libraries(${target} libmarqov)
set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )

endforeach(target)

if(MPI_CXX_FOUND)
if(CMAKE_VERSION VERSION_LESS 3.10)
    target_include_directories(mpiMARQOV PUBLIC ${MPI_INCLUDE_PATH})
    target_link_libraries(mpiMARQOV ${MPI_C_LIBRARIES})
else()
    target_include_directories(mpiMARQOV PUBLIC ${MPI_CXX_INCLUDE_DIRS})
    target_link_libraries(mpiMARQOV ${MPI_CXX_LIBRARIES})
    set_target_properties(mpiMARQOV PROPERTIES LINK_FLAGS " ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    target_compile_options(mpiMARQOV PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
    target_compile_definitions(mpiMARQOV PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS})

endif()
endif(MPI_CXX_FOUND)

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h)

# Let's check whether we have SWIG and can hence build the language bindings.

if(CMAKE_VERSION VERSION_GREATER 3.7)
    find_package(SWIG 3.0 OPTIONAL_COMPONENTS python)

    if(SWIG_FOUND)
      cmake_policy(SET CMP0086 NEW)
      cmake_policy(SET CMP0078 NEW)
      FIND_PACKAGE(PythonLibs)
      INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
      INCLUDE(${SWIG_USE_FILE})
      message("SWIG found: ${SWIG_EXECUTABLE}.")
      set_property(SOURCE libshowcase.i PROPERTY CPLUSPLUS ON)

      message("PYTHONLIBS_VERSION_STRING: ${PYTHONLIBS_VERSION_STRING}")
      message("CMAKE_SWIG_FLAGS: ${CMAKE_SWIG_FLAGS}")

      swig_add_library(showcase TYPE SHARED LANGUAGE python SOURCES libshowcase.i)
      SWIG_LINK_LIBRARIES(showcase ${PYTHON_LIBRARIES})
    endif()


    # specific language queries not supported in cmake 3-16...
    #  if(NOT SWIG_python_FOUND)
    #    message(WARNING "SWIG python bindings cannot be generated")
    #  endif()
endif()
