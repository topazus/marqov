find_package(HDF5 COMPONENTS CXX)
# find_package(OpenMP)
find_package(Threads)

if(NOT HDF5_FOUND)
include(ExternalProject)
    set(hdf5_local_major_version "1.12")
    set(hdf5_local_minor_version "0")
    find_package(ZLIB REQUIRED)
    #HDF5 also requires in its configure step a C++ compiler
    ExternalProject_Add(HDF5_LOCAL
        URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_local_major_version}/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}/src/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}.tar.gz
        PREFIX hdf5_local
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-build-mode=production --enable-cxx
        BUILD_COMMAND make
        INSTALL_COMMAND make install
    )
    set(HDF5_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/include)
    set(HDF5_Fortran_LIBRARIES ${PROJECT_BINARY_DIR}/hdf5_local/lib/libhdf5_fortran.so;${ZLIB_LIBRARIES})
    set(HDF5_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/lib/)
    set(HDF5_Fortran_COMPILER_EXECUTABLE ${PROJECT_BINARY_DIR}/hdf5_local/bin/h5fc)
    #Try again to find HDF5
    set(ENV{HDF5_ROOT} ${PROJECT_BINARY_DIR}/hdf5_local)
    find_package(HDF5 COMPONENTS CXX)
endif()


INCLUDE(CheckFunctionExists)

SET(MARQOVTARGETS MARQOV mpiMARQOV)

foreach(target ${MARQOVTARGETS})
add_executable (${target} registry.cpp marqov.cpp main.cpp)
target_include_directories(${target} PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(${target} ${HDF5_LIBRARIES} Threads::Threads)
set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
# target_compile_features(${target} PUBLIC cxx_return_type_deduction)



#add_definitions("-fopenmp")
#set_target_properties(MARQOV PROPERTIES LINK_FLAGS " ${LINK_FLAGS} -fopenmp")

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
#set_target_properties(mpiMARQOV PROPERTIES LINK_FLAGS " ${LINK_FLAGS} -fopenmp")
endif()
#SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif(MPI_CXX_FOUND)

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

