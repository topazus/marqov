find_package(HDF5 COMPONENTS CXX)
# find_package(OpenMP)
find_package(Threads)

set(MYHDF5INCLUDES ${HDF5_INCLUDE_DIRS})
set(MYHDF5LIBS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(NOT HDF5_FOUND)
    include(FetchContent)
    FetchContent_Declare(libhdf5
    GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
    GIT_TAG hdf5-1_12_0)
    FetchContent_GetProperties(libhdf5)
    message("teset1")
    if(NOT libhdf5_POPULATED)
	FetchContent_Populate(libhdf5)
	message("teset2")
        # Make subproject to use 'BUILD_SHARED_LIBS=ON' setting.
        set(HDF5_EXTERNALLY_CONFIGURED 1)
        set(HDF5_EXPORTED_TARGETS "hdf5_cpp-static")
        set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Build FORTRAN support")
        set(HDF5_BUILD_CPP_LIB ON CACHE BOOL "Build HDF5 C++ Library")
        set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "Build HIGH Level HDF5 Library")
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Shared Libraries")
        set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Build HDF5 Tools")
        set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Build HDF5 Library Examples")
        set(BUILD_TESTING OFF CACHE BOOL "Build HDF5 Unit Testing")
        set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "Enable parallel build (requires MPI)")

#        set (HDF5_EXTERNAL_LIB_PREFIX "marqovhdf5")
	message(STATUS "libhdf5 source dir: ${libhdf5_SOURCE_DIR}")
	message(STATUS "libhdf5 binary dir: ${libhdf5_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR} ${TARGET}")
        add_subdirectory(${libhdf5_SOURCE_DIR} ${libhdf5_BINARY_DIR})
        set(HDF5_DIR ${libhdf5_BINARY_DIR} CACHE PATH "libhdf5 path"  FORCE)
#        set(ENV{HDF5_ROOT} ${PROJECT_BINARY_DIR}/hdf5_local)
#        find_package(HDF5 COMPONENTS CXX NO_MODULE)
    endif()
#    find_package(HDF5 COMPONENTS CXX)
set(MYHDF5INCLUDES ${HDF5_INCLUDE_DIR} ${libhdf5_BINARY_DIR})
set(MYHDF5LIBS hdf5_cpp-static)

#include(ExternalProject)
#    set(hdf5_local_major_version "1.12")
#    set(hdf5_local_minor_version "0")
#    find_package(ZLIB)
#    #HDF5 also requires in its configure step a C++ compiler
#    ExternalProject_Add(libhdf5
#        URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${hdf5_local_major_version}/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}/src/hdf5-${hdf5_local_major_version}.${hdf5_local_minor_version}.tar.gz
#    PREFIX hdf5_local
#        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-cxx
#        BUILD_COMMAND make
#        INSTALL_COMMAND make install
#    URL_HASH SHA256=5f9a3ee85db4ea1d3b1fa9159352aebc2af72732fc2f58c96a3f0768dba0e9aa
#    )
#    target_include_directories(libhdf5
#)
#    set(HDF5_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/include)
#    set(HDF5_CXX_LIBRARIES ${PROJECT_BINARY_DIR}/hdf5_local/lib/libhdf5_cxx.so;${ZLIB_LIBRARIES})
#    set(HDF5_LIBRARY_DIRS ${PROJECT_BINARY_DIR}/hdf5_local/lib/)
#    set(HDF5_CXX_COMPILER_EXECUTABLE ${PROJECT_BINARY_DIR}/hdf5_local/bin/h5c++)
    #Try again to find HDF5
#    set(ENV{HDF5_ROOT} ${PROJECT_BINARY_DIR}/hdf5_local)
#    find_package(HDF5 COMPONENTS CXX)



endif()


INCLUDE(CheckFunctionExists)

SET(MARQOVTARGETS MARQOV mpiMARQOV)

foreach(target ${MARQOVTARGETS})
add_executable (${target} registry.cpp marqov.cpp main.cpp)
target_include_directories(${target} PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(${target} ${MYHDF5LIBS} Threads::Threads)
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

