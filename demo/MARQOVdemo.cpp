/* MARQOV - A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

// MARQOV
#include "../src/libmarqov/libmarqov.h"
#include "../src/libmarqov/util/startup.h"



using namespace MARQOV;

int main(int argc, char* argv[])
{
#ifdef MPIMARQOV
    int threadingsupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &threadingsupport);
	//FIXME: maybe we get by with one level less.
    if(threadingsupport < MPI_THREAD_SERIALIZED)
    {
        std::cout << "[MARQOV::main] Couldn't initialize MPI! ";
		std::cout << "Requested threading level not supported." << std::endl;
        return -1;
    }
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) 
	{
#endif
	print_startup_message();
#ifdef MPIMARQOV
    }
#endif
	// read registry
	// create config folder and select.ini file if not available in the working dir
	RegistryDB registry;
	check_registry_availability(registry, "Ising");
	// run the actual simulation
//	selectsim(registry);

#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
