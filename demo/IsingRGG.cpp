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

#include "../src/libmarqov/libmarqov.h"
#include "../src/libmarqov/util/startup.h"
#include "../src/libmarqov/util/registry.h"
#include "../src/lattice/random_geometric_graph.h"
#include "../src/hamiltonian/Ising.h"

using namespace MARQOV;

int main(int argc, char* argv[])
{
	// MPI startup
	// -----------
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


	// Registry
	// --------
	const auto ham = "IsingRGG";
	const auto configfile = "IsingRGG.ini";

	// We utilize a small helper library of ours, the registry which
	// reads Windows .ini style files to read config files
	RegistryDB registry;
	check_registry_availability(registry, ham);
	check_registry_file_exists(registry, ham);
    printInfoandcheckreplicaconfig(registry, ham);

    std::string outbasedir = registry.Get<std::string>(configfile, "IO", "outdir" );
	makeDir(outbasedir);

	// Parameters
	// ----------
	const auto nL            = registry.Get<std::vector<int>>(configfile, ham, "L" );
	const auto dim           = registry.Get<int>(configfile, ham, "dim" );
	const int nreplicas      = registry.Get<int>(configfile, ham, "rep" );
	const auto nclusteramp   = registry.Get<double>(configfile, "MC", "nclusteramp");
	const auto nclusterexp   = registry.Get<int>(configfile, "MC", "nclusterexp");
	const auto nmetro        = registry.Get<int>(configfile, "MC", "nmetro");
	const auto warmupsteps   = registry.Get<int>(configfile, "MC", "warmupsteps");
	const auto measuresteps  = registry.Get<int>(configfile, "MC", "measuresteps");
	const auto nthreads      = number_of_threads_per_node(registry, configfile);

	// Typedefs
	// --------
	typedef Ising<int> Hamiltonian;
	typedef RandomGeometricGraph<Poissonian> Lattice;
	typedef std::knuth_b RNG;


	// Hamltonian parameters
	// ---------------------
	auto beta = registry.Get<std::vector<double> >(configfile, ham, "beta");
	auto J    = registry.Get<std::vector<double> >(configfile, ham, "J");
	auto hp = cart_prod(beta, J);


	// Loop over lattice sizes
	for (std::size_t j=0; j<nL.size(); j++)
	{
		// Print status and prepare directory
		int L = nL[j];
		std::cout << std::endl << "L = " << L << std::endl << std::endl;
		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		makeDir(outpath);


		// Monte Carlo parameters
		// ----------------------
		MARQOV::Config mp(outpath);
		mp.setnmetro(nmetro); // number of Metropolis sweeps per EMCS
		mp.setncluster(int(nclusteramp*pow(L,nclusterexp))); // number of Wolff updates per EMCS
		mp.setwarmupsteps(warmupsteps); // number of EMCS for warmup
		mp.setgameloopsteps(measuresteps); // number of EMCS for production


		// Geometry parameters
		// -------------------
		// search radius for RGG
		double no_nn = 6;
		double search_radius = std::cbrt((3.0/M_PI*no_nn*1.0/4.0))/double(L); // only for 3D


		// Schedule and Run
		// ----------------

		// Bundle the lattice, the MC parameters and the Hamiltonian parameters
		auto paramsets  = finalize_parameter(std::make_tuple(L, dim, search_radius), mp, hp);

		// Duplicate everything for the amount of replicas
		auto rparamsets = replicator(paramsets, nreplicas);

		// Instantiate the scheduler which waits for new threads
		// (makeScheduler can figure out a lot from one set of parameters)
		auto sched = makeScheduler<Hamiltonian, Lattice> (rparamsets[0], 1, nthreads);

		// Submit parameter sets to the scheduler 
		for (auto p: rparamsets) sched.createSimfromParameter(p, defaultfilter);

		// Run!
		sched.start();
	}

#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
