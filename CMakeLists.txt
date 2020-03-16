cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(marqov CXX)

add_subdirectory(src)

find_package(MPI REQUIRED CXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
target_link_libraries(marqov       ${MPI_C_LIBRARIES})
SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

#add_subdirectory(testsuite)
#enable_testing()
