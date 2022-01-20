/* MARQOV - A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2022, The MARQOV Project
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
#include <tuple>
#include <cstdlib>

#include "../src/libmarqov/libmarqov.h"
#include "../src/libmarqov/util/startup.h"
#include "../src/libmarqov/util/registry.h"
#include "../src/libmarqov/util/regularlatticeloop.h"
#include "../src/hamiltonian/Phi4.h"

using namespace MARQOV;

int main(int argc, char* argv[])
{
	// MPI startup
	#ifdef MPIMARQOV
    int threadingsupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &threadingsupport);
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

	// -----------------------------------------

	const auto ham = "Phi4";
	const auto configfile = "Phi4.ini";

	// Load config
	RegistryDB registry;
	check_registry_availability(registry, ham);
	check_registry_file_exists(registry, ham);
    printInfoandcheckreplicaconfig(registry, ham);
	
	// Prepare output folder
	auto outbasedir = registry.Get<std::string>(configfile, "IO", "outdir" );
	makeDir(outbasedir);

	// Parameters
	auto beta   = registry.Get<std::vector<double> >(configfile, ham, "beta");
	auto lambda = registry.Get<std::vector<double> >(configfile, ham, "lambda");
	auto mass   = registry.Get<std::vector<double> >(configfile, ham, "mass");

	// we need "beta" as an explicit parameter in the Hamiltonian
	// this requires some gymnastics when using "cart_prod"
	std::vector<double> dummy = {0.0};
	auto parameters = cart_prod(beta, dummy, lambda, mass);
	for (std::size_t i=0; i<parameters.size(); i++)
		std::get<1>(parameters[i]) = std::get<0>(parameters[i]);


	// Execute the actual simulations
	RegularLatticeLoop<Phi4<double,double>>(registry, outbasedir, parameters, defaultfilter);


#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
