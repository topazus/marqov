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

// MARQOV
#include "../src/libmarqov/libmarqov.h"
#include "../src/libmarqov/util/startup.h"

// Lattices
#include "../src/lattice/random_geometric_graph.h"

// Hamiltonians
#include "../src/hamiltonian/Ising.h"

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


	// Input 
	// -----
	const auto ham = "IsingRGG";
	const auto configfile = "IsingRGG.ini";
	// We utilize a small helper library of ours, the registry which
	// reads Windows .ini style files to read config files
	RegistryDB registry;
	check_registry_availability(registry, ham);
	check_registry_file_exists(registry, ham);


	// Parameters
    std::string outbasedir = registry.Get<std::string>(configfile, "IO", "outdir" );
    tidyupoldsims(outbasedir);

	auto nreplicas  = registry.Get<std::vector<int>>(configfile, ham, "rep" );
	const auto nL   = registry.Get<std::vector<int>>(configfile, ham, "L" );
	const auto dim  = registry.Get<int>(configfile, ham, "dim" );
    printInfoandcheckreplicaconfig(registry, ham);


	// Number of threads
	int nthreads = 0;
	try 
	{
		nthreads = registry.template Get<int>(configfile, "General", "threads_per_node" );
	}
	catch (const Registry_Key_not_found_Exception&) 
	{
		std::cout<<"threads_per_node not set -> automatic"<<std::endl;
	}


	// Replicas
	if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

	// Physical parameters
	auto beta = registry.Get<std::vector<double> >(configfile, ham, "beta");
	auto J    = registry.Get<std::vector<double> >(configfile, ham, "J");
	auto hp = cart_prod(beta, J);

	// Typedefs
	typedef Ising<int> Hamiltonian;
	typedef RandomGeometricGraph<Poissonian> Lattice;
	typedef std::knuth_b RNG;

    typedef std::tuple<std::tuple<int, int, double>, MARQOV::Config, typename decltype(hp)::value_type > ParameterType;
	typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType, RNG>::MarqovScheduler SchedulerType;


	// Lattice size loop
	for (std::size_t j=0; j<nL.size(); j++)
	{
		// init scheduler
		SchedulerType sched(1, nthreads);

		// prepare output
		int L = nL[j];
		std::cout << std::endl << "L = " << L << std::endl << std::endl;
		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		makeDir(outpath);

		// Monte Carlo parameters
		MARQOV::Config mp(outpath);
		mp.setnmetro(5);
		mp.setncluster(15);
		mp.setwarmupsteps(150);
		mp.setgameloopsteps(300);

		// search radius for RGG
		double no_nn = 6;
		double search_radius = std::cbrt((3.0/M_PI*no_nn*1.0/4.0))/double(L); // only for 3D

		// form parameter triple with lattice parameters and replicate
		auto params  = finalize_parameter(std::make_tuple(L, dim, search_radius), mp, hp);
		auto rparams = replicator(params, nreplicas[j]);

		// feed scheduler
		for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);

		// run!
		sched.start();
	}

#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
