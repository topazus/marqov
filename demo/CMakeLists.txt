INCLUDE(CheckFunctionExists)

#add_subdirectory(../src/libmarqov)

# configuration for the simple ising test case
add_executable(mysimpleising mysimpleising.cpp)
target_include_directories(mysimpleising PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(mysimpleising libMARQOV)
set_target_properties(mysimpleising PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

# configuration for the Heisenberg test case with threading
add_executable(mysimpleheisenbergthreaded mysimpleheisenbergthreaded.cpp)
target_include_directories(mysimpleheisenbergthreaded PUBLIC ${MYHDF5INCLUDES})
target_link_libraries(mysimpleheisenbergthreaded libMARQOV)
set_target_properties(mysimpleheisenbergthreaded PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)


# configuration of the main parts.

SET(MARQOVTARGETS MARQOVdemo mpiMARQOVdemo)

foreach(target ${MARQOVTARGETS})
    add_executable (${target} ../src/libmarqov/util/registry.cpp MARQOVdemo.cpp)
    target_link_libraries(${target} libMARQOV)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )

endforeach(target)

find_package(PCL 1.2 COMPONENTS common kdtree REQUIRED)
link_directories(${PCL_LIBRARY_DIRS}) 
add_definitions(${PCL_DEFINITIONS})
include_directories(${PCL_INCLUDE_DIRS})
target_link_libraries(MARQOVdemo           ${PCL_LIBRARIES})




find_package(MPI)
if(MPI_CXX_FOUND)
    target_compile_options(mpiMARQOVdemo PUBLIC -DMPIMARQOV)
    if(CMAKE_VERSION VERSION_LESS 3.10)
	target_include_directories(mpiMARQOVdemo PUBLIC ${MPI_INCLUDE_PATH})
	target_link_libraries(mpiMARQOVdemo ${MPI_CXX_LIBRARIES})
	target_compile_options(mpiMARQOVdemo PUBLIC ${MPI_CXX_COMPILE_FLAGS})
    else()
	target_include_directories(mpiMARQOVdemo PUBLIC ${MPI_CXX_INCLUDE_DIRS})
	target_link_libraries(mpiMARQOVdemo ${MPI_CXX_LIBRARIES})
	set_target_properties(mpiMARQOVdemo PROPERTIES LINK_FLAGS " ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
	target_compile_options(mpiMARQOVdemo PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
	target_compile_definitions(mpiMARQOVdemo PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS})

    endif()
endif(MPI_CXX_FOUND)

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)

CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)

CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/libmarqov/util/config.h)
