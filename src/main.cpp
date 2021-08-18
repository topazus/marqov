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
#include "libmarqov/libmarqov.h"

// Lattices
#include "lattice/regular_hypercubic.h"
#include "lattice/constant_coordination.h"
#include "lattice/regular_random_bond.h"
#include "lattice/simple_bipartite.h"


// Hamiltonians
#include "hamiltonian/Heisenberg.h"
#include "hamiltonian/Ising.h"
#include "hamiltonian/Potts.h"
#include "hamiltonian/Phi4.h"
#include "hamiltonian/BlumeCapel.h"
#include "hamiltonian/BlumeEmeryGriffiths.h"
#include "hamiltonian/XXZAntiferro.h"
#include "hamiltonian/XXZAntiferroSingleAniso.h"
#include "hamiltonian/AshkinTeller.h"
#include "hamiltonian/EdwardsAndersonIsing.h"
//#include "hamiltonian/Ssh.h" // seperate branch
#include "hamiltonian/BlumeCapelBipartite.h"
#include "hamiltonian/AshkinTeller.h"

using namespace MARQOV;

bool startswith(const std::string longword, const std::string shortword)
{
	if (longword.rfind(shortword, 0) == 0) return true;
	else return false; 
}


std::string selectsim_startup(RegistryDB& registry)
{
	const auto ham        = registry.Get<std::string>("select.ini", "General", "Hamiltonian" );
	const auto dim 	      = registry.Get<int>(ham+".ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
	const auto nreplicass = registry.Get<std::string>(ham+".ini", ham, "rep" );
	const auto nL  	      = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
	const auto nLs 	      = registry.Get<std::string>(ham+".ini", ham, "L" );

	cout << endl;
	cout << "Hamiltonian: \t" << ham << endl;
	cout << "Dimension: \t" << dim << endl;
	cout << "Lattice sizes:\t" << nLs << endl;
	cout << "Replicas:\t" << nreplicass << endl;

	if ((nreplicas.size() != nL.size()) && (nreplicas.size() != 1)) throw std::invalid_argument("invalid replica configuration!");


	return ham;
}





// ---------------------------------------

void selectsim()
{

	RegistryDB registry("../src/config", "ini");
	const auto ham = selectsim_startup(registry);

	std::string outbasedir = registry.Get<std::string>(ham+".ini", "IO", "outdir" );

	// delete previous output // fixme: don't do that by default!
	std::string command;
	command = "rm -r " + outbasedir;
	system(command.c_str());
	makeDir(outbasedir);


	// ----------------- select simulation ------------------

	if (startswith(ham,"Ising"))
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto parameters = cart_prod(beta, J);

 		RegularLatticeLoop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
	}

	else if (startswith(ham,"Potts"))
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto parameters = cart_prod(beta, J);

		switch(registry.Get<int>(ham+".ini", ham, "q"))
		{
			case 3:
				RegularLatticeLoop<Potts<3>>(registry, outbasedir, parameters, defaultfilter);
				break;
			case 4:
				RegularLatticeLoop<Potts<4>>(registry, outbasedir, parameters, defaultfilter);
				break;
			case 6:
				RegularLatticeLoop<Potts<6>>(registry, outbasedir, parameters, defaultfilter);
				break;
			case 8:
				RegularLatticeLoop<Potts<8>>(registry, outbasedir, parameters, defaultfilter);
				break;
			default:
				std::cout<<"[MARQOV::main] Potts: unsupported q!";
		}
	}

	else if (ham == "AshkinTeller")
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto K    = registry.Get<std::vector<double> >(ham+".ini", ham, "K");
		auto parameters = cart_prod(beta, J, K);

		RegularLatticeLoop<AshkinTeller>(registry, outbasedir, parameters, defaultfilter);
	}


	else if (ham == "Heisenberg")
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto parameters = cart_prod(beta, J);

		RegularLatticeLoop<Heisenberg<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}



	else if (ham == "Phi4")
	{

		auto beta   = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto lambda = registry.Get<std::vector<double> >(ham+".ini", ham, "lambda");
		auto mass   = registry.Get<std::vector<double> >(ham+".ini", ham, "mass");
		
		// we need "beta" as an explicit parameter in the Hamiltonian
		// this requires some gymnastics ...
		std::vector<double> dummy = {0.0};
		auto parameters = cart_prod(beta, dummy, lambda, mass);
		for (std::size_t i=0; i<parameters.size(); i++) 
			std::get<1>(parameters[i]) = std::get<0>(parameters[i]);
		
		RegularLatticeLoop<Phi4<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}



	else if (ham == "BlumeCapel")
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto D    = registry.Get<std::vector<double> >(ham+".ini", ham, "D");
		auto parameters = cart_prod(beta, J, D);
		
		RegularLatticeLoop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
	}



	else if (ham == "BlumeEmeryGriffiths")
	{
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto D    = registry.Get<std::vector<double> >(ham+".ini", ham, "D");
		auto K    = registry.Get<std::vector<double> >(ham+".ini", ham, "K");
		auto parameters = cart_prod(beta, J, D, K);
		
		RegularLatticeLoop<BlumeEmeryGriffiths<int>>(registry, outbasedir, parameters, defaultfilter);
	}



	else if (ham == "XXZAntiferro")
	{
		auto beta     = registry.Get<std::vector<double>>(ham+".ini", ham, "beta");
		auto extfield = registry.Get<std::vector<double>>(ham+".ini", ham, "extfield");
		auto aniso    = registry.Get<std::vector<double>>(ham+".ini", ham, "aniso");
		auto parameters = cart_prod(beta, aniso, extfield);
		
		RegularLatticeLoop<XXZAntiferro<double>>(registry, outbasedir, parameters, defaultfilter);
	}



	else if (ham == "XXZAntiferroSingleAniso")
	{
		auto beta        = registry.Get<std::vector<double>>(ham+".ini", ham, "beta");
		auto extfield    = registry.Get<std::vector<double>>(ham+".ini", ham, "extfield");
		auto aniso       = registry.Get<std::vector<double>>(ham+".ini", ham, "aniso");
		auto singleaniso = registry.Get<std::vector<double>>(ham+".ini", ham, "singleaniso");
		auto parameters = cart_prod(beta, extfield, aniso, singleaniso);

		RegularLatticeLoop<XXZAntiferroSingleAniso<double>>(registry, outbasedir, parameters, xxzfilter);
	}


	else if (startswith(ham, "EdwardsAnderson-Ising"))
	{
		// Parameters
		auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
		const auto nL   = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
		const auto dim  = registry.Get<int>(ham+".ini", ham, "dim" );
	
	
		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>(ham+".ini", ham, "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}

		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		// Physical parameters
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto hp = cart_prod(beta, J);
        

		// Typedefs
		typedef EdwardsAndersonIsing<int> Hamiltonian;
		typedef RegularRandomBond<GaussianPDF> Lattice;
        //typedef RegularRandomBond<BimodalPDF> Lattice;

        typedef std::tuple<std::tuple<int, int>, MARQOV::Config, typename decltype(hp)::value_type > ParameterType;
		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType, std::ranlux48>::MarqovScheduler SchedulerType;

		SchedulerType sched(1, nthreads);

		// Lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
			MARQOV::Config mp(outpath);
			mp.setnmetro(50);
			mp.setncluster(0);
			mp.setwarmupsteps(200);
			mp.setgameloopsteps(1000);

			// lattice parameters
// 			auto lp = std::make_tuple(L, dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter(std::make_tuple(L, dim), mp, hp);//this particular form is required to happify PGI-19.10
            auto rparams = replicator(params, nreplicas[j]);

			// schedule simulations
 			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);
		}
 		sched.start(); // run!
	}

	else if (ham == "IsingCC")
	{
		// Parameters
		auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
		const auto nL   = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
		const auto dim  = registry.Get<int>(ham+".ini", ham, "dim" );
	
	
		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>(ham+".ini", ham, "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}


		// Replicas
		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		// Physical parameters
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto hp = cart_prod(beta, J);

		// Typedefs
		typedef Ising<int> Hamiltonian;
		typedef ConstantCoordinationLattice<Poissonian> Lattice;

		// Lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
			MARQOV::Config mp(outpath);
			mp.setnmetro(5);
			mp.setncluster(15);
			mp.setwarmupsteps(500);
			mp.setgameloopsteps(1500);

			// form parameter triple with lattice parameters and replicate
			auto params  = finalize_parameter(std::make_tuple(L, dim), mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

			// init scheduler            
            auto sched = makeScheduler<Hamiltonian, Lattice, std::knuth_b>(params[0], 1, nthreads);
			// feed scheduler
			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);

			// run!
			sched.start();
		}
	}

	else if (ham == "BlumeCapelBipartite")
	{
		// Parameters
		auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
		const auto nL   = registry.Get<std::vector<int>>(ham+".ini", ham,  "L" );
		const auto dim  = registry.Get<int>(ham+".ini", ham, "dim" );
	
	
		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>(ham+".ini", "ham", "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}


		// Replicas
		if (nreplicas.size() == 1) { for (std::size_t i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }


		// import parameters
		auto beta = registry.Get<std::vector<double> >(ham+".ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >(ham+".ini", ham, "J");
		auto DA   = registry.Get<std::vector<double> >(ham+".ini", ham, "DA");
		auto DB   = registry.Get<std::vector<double> >(ham+".ini", ham, "DB");
		auto hp = cart_prod(beta, J, DA, DB);

		typedef BlumeCapelBipartite<int> Hamiltonian;
		typedef SimpleBipartite Lattice;

        typedef typename std::tuple<SimpleBipartite&, MARQOV::Config, std::tuple<double, double, double, double> > ParameterType;
		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;

		// Prepare Geometry
		std::vector<SimpleBipartite> latts;
		for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j], dim);
	    
	
		// Init Scheduler
		SchedulerType sched(1, nthreads);
	    
	
		// Lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
	
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
	
			MARQOV::Config mp(outpath);
			mp.setnmetro(10);
			mp.setncluster(int(L/2));
			mp.setwarmupsteps(200);
			mp.setgameloopsteps(500);
			
			makeDir(mp.outpath);
			
			// set up and execute        
			Lattice& latt = latts[j];
			auto params = finalize_parameter(latt, mp, hp);
            auto rparams = replicator(params, nreplicas[j]);
	
			// feed the scheduler
			for(auto p: rparams) sched.createSimfromParameter(p, defaultfilter);
		}
		sched.start();
	}
}

int main(int argc, char* argv[])
{

    std::cout<<"MARQOV Copyright (C) 2020-2021, The MARQOV Project contributors"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;

	// read config files
	RegistryDB registry("../src/config", "ini");


#ifdef MPIMARQOV
    int threadingsupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &threadingsupport);//FIXME: maybe we get by with one level less.
    if(threadingsupport < MPI_THREAD_SERIALIZED)
    {
        std::cout<<"couldn't initialize MPI! threading level not supported."<<std::endl;
        return -1;
    }

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
    }
#endif

    selectsim();

#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
