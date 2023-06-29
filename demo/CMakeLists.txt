INCLUDE(CheckFunctionExists)

add_executable(IsingSimple IsingSimple.cpp)
add_executable(HeisenbergSimpleThreaded HeisenbergSimpleThreaded.cpp)

SET(DEMOS Ising Heisenberg Phi4 Potts BlumeCapel XXZAntiferroSingleAniso BlumeCapelBipartiteLattice AshkinTellerThreeColor BlumeEmeryGriffiths)
foreach(demo ${DEMOS})
add_executable(${demo} ${demo}.cpp ../src/libmarqov/util/registry.cpp)
install(TARGETS ${demo}
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/marqov/demo)

add_executable(${demo}.mpi ${demo}.cpp ../src/libmarqov/util/registry.cpp)
install(TARGETS ${demo}.mpi
          DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/marqov/demo)
endforeach(demo)

# lists make things easier ...
SET(DEMOTARGETS IsingSimple HeisenbergSimpleThreaded ${DEMOS})

# SET(DEMOTARGETSMPI ${DEMOS})
foreach(demo ${DEMOS})
SET(DEMOTARGETSMPI "${DEMOTARGETSMPI};${demo}.mpi")
endforeach()
# list(TRANSFORM DEMOTARGETSMPI APPEND ".mpi")

# compile IsingRGG demo only if PCL library is installed
find_package(PCL 1.2 COMPONENTS common kdtree)
if(PCL_FOUND)
	add_executable(IsingRGG IsingRGG.cpp ../src/libmarqov/util/registry.cpp)
	add_executable(IsingRGG.mpi IsingRGG.cpp ../src/libmarqov/util/registry.cpp)
	SET(DEMOTARGETS "${DEMOTARGETS};IsingRGG")
	SET(DEMOTARGETSMPI "${DEMOTARGETSMPI};IsingRGG.mpi")

	target_link_libraries(IsingRGG ${PCL_LIBRARIES})
        target_compile_definitions(IsingRGG PUBLIC ${PCL_DEFINITIONS})
        target_include_directories(IsingRGG PUBLIC ${PCL_INCLUDE_DIRS})
	target_link_libraries(IsingRGG.mpi ${PCL_LIBRARIES})
        target_compile_definitions(IsingRGG.mpi PUBLIC ${PCL_DEFINITIONS})
        target_include_directories(IsingRGG.mpi PUBLIC ${PCL_INCLUDE_DIRS})
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

# Set common options for all targets
SET(ALLDEMOS "${DEMOTARGETS};${DEMOTARGETSMPI}")
foreach(target ${ALLDEMOS})
    target_link_libraries(${target} libMARQOV)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)
	target_include_directories(${target} PUBLIC ${MYHDF5INCLUDES})
    set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMPILE_FLAGS " ${COMPILE_FLAGS} -DHAVE_CONFIG_H -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_BRANCH=\"${GIT_BRANCH}\"" )
endforeach(target)

# groups
# this allows to do "make demo" and "make demo.mpi"
add_custom_target(demo DEPENDS ${DEMOTARGETS})
add_custom_target(demo.mpi DEPENDS ${DEMOTARGETSMPI})

CHECK_FUNCTION_EXISTS(scandir HAVE_SCANDIR)
CHECK_FUNCTION_EXISTS(alphasort HAVE_ALPHASORT)
CHECK_FUNCTION_EXISTS(readdir_r HAVE_READDIR_R)
