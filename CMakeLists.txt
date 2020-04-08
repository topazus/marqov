cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(marqov CXX)

if(CMAKE_VERSION VERSION_LESS "3.10.0")
find_package(MPI)
else()
find_package(MPI COMPONENTS CXX)

endif()

add_subdirectory(src)

#add_subdirectory(testsuite)
#enable_testing()
