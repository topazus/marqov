cmake_minimum_required(VERSION 3.0.0)
Project(MARQOV.tests CXX)

find_package (Eigen3 QUIET)

add_subdirectory(components)
# add_subdirectory(systems)

if (TARGET Eigen3::Eigen)
    add_executable(4-eigentest 4-eigentest.cpp)

    target_link_libraries(4-eigentest libMARQOV)
    target_link_libraries(4-eigentest Eigen3::Eigen)
    target_include_directories(4-eigentest PUBLIC ${CMAKE_SOURCE_DIR}/src/libmarqov)
endif (TARGET Eigen3::Eigen)


enable_testing()
if (TARGET Eigen3::Eigen)
    add_test(4-eigentest 4-eigentest)
endif (TARGET Eigen3::Eigen)