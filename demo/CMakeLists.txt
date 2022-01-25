INCLUDE(CheckFunctionExists)

add_executable(IsingSimple IsingSimple.cpp)

add_executable(HeisenbergSimpleThreaded HeisenbergSimpleThreaded.cpp)

add_executable(Ising Ising.cpp ../src/libmarqov/util/registry.cpp)
add_executable(Ising.mpi Ising.cpp ../src/libmarqov/util/registry.cpp)

add_executable(Heisenberg Heisenberg.cpp ../src/libmarqov/util/registry.cpp)
add_executable(Heisenberg.mpi Heisenberg.cpp ../src/libmarqov/util/registry.cpp)

add_executable(Phi4 Phi4.cpp ../src/libmarqov/util/registry.cpp)
add_executable(Phi4.mpi Phi4.cpp ../src/libmarqov/util/registry.cpp)

add_executable(BlumeCapel BlumeCapel.cpp ../src/libmarqov/util/registry.cpp)
add_executable(BlumeCapel.mpi BlumeCapel.cpp ../src/libmarqov/util/registry.cpp)

add_executable(XXZAntiferroSingleAniso XXZAntiferroSingleAniso.cpp ../src/libmarqov/util/registry.cpp)
add_executable(XXZAntiferroSingleAniso.mpi XXZAntiferroSingleAniso.cpp ../src/libmarqov/util/registry.cpp)

add_executable(BlumeCapelBipartiteLattice BlumeCapelBipartiteLattice.cpp ../src/libmarqov/util/registry.cpp)
add_executable(BlumeCapelBipartiteLattice.mpi BlumeCapelBipartiteLattice.cpp ../src/libmarqov/util/registry.cpp)

add_executable(AshkinTellerThreeColor AshkinTellerThreeColor.cpp ../src/libmarqov/util/registry.cpp)
add_executable(AshkinTellerThreeColor.mpi AshkinTellerThreeColor.cpp ../src/libmarqov/util/registry.cpp)

add_executable(BlumeEmeryGriffiths BlumeEmeryGriffiths.cpp ../src/libmarqov/util/registry.cpp)
add_executable(BlumeEmeryGriffiths.mpi BlumeEmeryGriffiths.cpp ../src/libmarqov/util/registry.cpp)


SET(DEMOTARGETS IsingSimple HeisenbergSimpleThreaded Ising Heisenberg Phi4 BlumeCapel XXZAntiferroSingleAniso BlumeCapelBipartiteLattice AshkinTellerThreeColor BlumeEmeryGriffiths)

SET(DEMOTARGETSMPI Ising.mpi Heisenberg.mpi Phi4.mpi BlumeCapel.mpi XXZAntiferroSingleAniso.mpi BlumeCapelBipartiteLattice.mpi AshkinTellerThreeColor.mpi BlumeEmeryGriffiths.mpi)


# link non-MPI targets
foreach(target ${DEMOTARGETS})
    target_link_libraries(${target} libMARQOV)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
	target_include_directories(${target} PUBLIC ${MYHDF5INCLUDES})
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
endforeach(target)

# link MPI targets
foreach(target ${DEMOTARGETSMPI})
    target_link_libraries(${target} libMARQOV)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
	target_include_directories(${target} PUBLIC ${MYHDF5INCLUDES})
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
endforeach(target)


# compile IsingRGG demo only if PCL library is installed
find_package(PCL 1.2 COMPONENTS common kdtree)
if(PCL_FOUND)
	add_executable(IsingRGG IsingRGG.cpp ../src/libmarqov/util/registry.cpp)
	add_executable(IsingRGG.mpi IsingRGG.cpp ../src/libmarqov/util/registry.cpp)

	link_directories(${PCL_LIBRARY_DIRS}) 
	add_definitions(${PCL_DEFINITIONS})
	include_directories(${PCL_INCLUDE_DIRS})

	target_link_libraries(IsingRGG ${PCL_LIBRARIES})
	target_link_libraries(IsingRGG.mpi ${PCL_LIBRARIES})

    target_link_libraries(IsingRGG libMARQOV)
    target_link_libraries(IsingRGG.mpi libMARQOV)

    set_target_properties(IsingRGG PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
	set_target_properties(IsingRGG.mpi PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

	target_include_directories(IsingRGG PUBLIC ${MYHDF5INCLUDES})
	target_include_directories(IsingRGG.mpi PUBLIC ${MYHDF5INCLUDES})

    set_target_properties(IsingRGG PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
    set_target_properties(IsingRGG.mpi PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
endif(PCL_FOUND)


# some models are MPI parallelized
find_package(MPI)
if(MPI_CXX_FOUND)
	foreach(target ${DEMOTARGETSMPI})
	    target_compile_options(${target} PUBLIC -DMPIMARQOV)
    	if(CMAKE_VERSION VERSION_LESS 3.10)
			target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
			target_link_libraries(${target} ${MPI_CXX_LIBRARIES})
			target_compile_options(${target} PUBLIC ${MPI_CXX_COMPILE_FLAGS})
    	else()
			target_include_directories(${target} PUBLIC ${MPI_CXX_INCLUDE_DIRS})
			target_link_libraries(${target} ${MPI_CXX_LIBRARIES})
			set_target_properties(${target} PROPERTIES LINK_FLAGS " ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
			target_compile_options(${target} PUBLIC ${MPI_CXX_COMPILE_OPTIONS})
			target_compile_definitions(${target} PUBLIC ${MPI_CXX_COMPILE_DEFINITIONS})
    	endif()
	endforeach(target)
endif(MPI_CXX_FOUND)

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)
CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)
CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)
