find_package(HDF5 COMPONENTS CXX)
add_executable (marqov main.cpp)
target_include_directories(marqov PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(marqov ${HDF5_LIBRARIES})
target_compile_features(marqov PUBLIC cxx_lambdas cxx_auto_type cxx_constexpr cxx_trailing_return_types cxx_return_type_deduction)

