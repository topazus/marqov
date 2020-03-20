find_package(HDF5 COMPONENTS CXX)
INCLUDE(CheckFunctionExists)
add_executable (marqov registry.cpp main.cpp)
target_include_directories(marqov PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(marqov ${HDF5_LIBRARIES})
target_compile_features(marqov PUBLIC cxx_lambdas cxx_auto_type cxx_constexpr cxx_trailing_return_types cxx_return_type_deduction)

add_definitions(-O3)

set_target_properties(marqov PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H" )

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

