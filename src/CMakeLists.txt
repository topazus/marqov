find_package(HDF5 COMPONENTS CXX)

find_package(Threads)

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
find_package(MPI)
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

