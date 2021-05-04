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

#include "timetracker.h"
#include "helpers.h"
#include "vectorhelpers.h"
#include "cartprod.h"
#include "registry.h"
#include "systemtools.h"
#include "replicate.h"
#include "svmath.h"
#include "filters.h"
#include "marqovscheduler.h"

// Geometry
#include "geometry/regular_lattice.h"
#include "geometry/ssh_lattice.h"
#include "geometry/grid.h"
#include "geometry/neighbourclass.h"
#include "geometry/io.h"


// Hamiltonians
#include "hamiltonian/Heisenberg.h"
#include "hamiltonian/Ising.h"
#include "hamiltonian/Phi4.h"
#include "hamiltonian/BlumeCapel.h"
#include "hamiltonian/XXZAntiferro.h"
#include "hamiltonian/XXZAntiferroSingleAniso.h"
#include "hamiltonian/AshkinTeller.h"
#include "hamiltonian/EdwardsAndersonIsing.h"
//#include "hamiltonian/Ssh.h" // seperate branch
#include "hamiltonian/BlumeCapelBipartite.h"

using namespace MARQOV;

// ---------------------------------------

template <class Hamiltonian, class Lattice, class Parameters, class Callable>
void Loop(const std::vector<Parameters>& params, Callable filter)
{
	typename GetSchedulerType<Hamiltonian, Lattice, Parameters >::MarqovScheduler sched(1);
	    
	for(auto p : params)
	sched.createSimfromParameter(p, filter);
	sched.start();
}

// ---------------------------------------


template <class Hamiltonian, class Params, class Callable>
void RegularLatticeLoop(RegistryDB& reg, const std::string outbasedir, const std::vector<Params>& hp, Callable filter)
{
	
	// Typedefs
	typedef decltype(finalize_parameter_pair(std::declval<MARQOV::Config>(), hp)) ParameterPairType;
	typedef typename ParameterPairType::value_type ParameterType;
	typedef typename GetSchedulerType<Hamiltonian, RegularHypercubic, ParameterType>::MarqovScheduler SchedulerType;


	// Parameters
	const auto name = reg.Get<std::string>("mc.ini", "General", "Hamiltonian" );
	auto nreplicas  = reg.Get<std::vector<int>>("mc.ini", name, "rep" );
	const auto nL   = reg.Get<std::vector<int>>("mc.ini", name, "L" );
	const auto dim  = reg.Get<int>("mc.ini", name, "dim" );


	// Number of threads
	int nthreads = 0;
	try 
	{
		nthreads = reg.template Get<int>("mc.ini", "General", "threads_per_node" );
	}
	catch (const Registry_Key_not_found_Exception&) 
	{
		std::cout<<"threads_per_node not set -> automatic"<<std::endl;
	}


	// Replicas
	if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }


	// Prepare Geometry
	std::vector<RegularHypercubic> latts;
	for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j], dim);
    

	// Init Scheduler
	SchedulerType sched(1, nthreads);
    


	for (std::size_t j=0; j<nL.size(); j++)
	{
		// prepare
		int L = nL[j];
		cout << endl << "L = " << L << endl << endl;

		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";

		MARQOV::Config mp(outpath);
		mp.setnsweeps(10);
		mp.setncluster(int(L/2));
		mp.setwarmupsteps(200);
		mp.setgameloopsteps(500);
		
		makeDir(mp.outpath);
		
		auto params = finalize_parameter_pair(mp, hp);
		auto rparams = replicator_pair(params, nreplicas[j]);
		
		// set up and execute        
		RegularHypercubic& latt = latts[j];
		auto f = [&filter, &latt, &outbasedir, L](auto p){return filter(latt, p);}; //partially apply filter

		// feed the scheduler
		for(auto p: rparams) sched.createSimfromParameter(p, f);
	}
	sched.start();
}


// ---------------------------------------

std::string selectsim_startup(RegistryDB& registry)
{
	const auto ham        = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
	const auto dim 	    = registry.Get<int>("mc.ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", ham, "rep" );
	const auto nreplicass = registry.Get<std::string>("mc.ini", ham, "rep" );
	const auto nL  	    = registry.Get<std::vector<int>>("mc.ini", ham, "L" );
	const auto nLs 	    = registry.Get<std::string>("mc.ini", ham, "L" );

	cout << endl;
	cout << "Hamiltonian: \t" << ham << endl;
	cout << "Dimension: \t" << dim << endl;
	cout << "Lattice sizes:\t" << nLs << endl;
	cout << "Replicas:\t" << nreplicass << endl;

	if ((nreplicas.size() != nL.size()) && (nreplicas.size() != 1)) throw std::invalid_argument("invalid replica configuration!");

	return ham;
}

// ---------------------------------------


void selectsim(RegistryDB& registry, std::string outbasedir, std::string logbasedir)
{

	auto ham = selectsim_startup(registry);



	// ----------------- select simulation ------------------



	if (startswith(ham, "Ising"))
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);
 		RegularLatticeLoop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "Heisenberg")
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);
		RegularLatticeLoop<Heisenberg<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "Phi4")
	{
		auto beta   = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto lambda = registry.Get<std::vector<double> >("mc.ini", ham, "lambda");
		auto mass   = registry.Get<std::vector<double> >("mc.ini", ham, "mass");
		
		// we need "beta" as an explicit parameter in the Hamiltonian
		// this requires some gymnastics ...
		std::vector<double> dummy = {0.0};
		auto parameters = cart_prod(beta, dummy, lambda, mass);
		for (std::size_t i=0; i<parameters.size(); i++) std::get<1>(parameters[i]) = std::get<0>(parameters[i]);
		
		write_logfile(registry, beta);
		RegularLatticeLoop<Phi4<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "BlumeCapel")
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto D    = registry.Get<std::vector<double> >("mc.ini", ham, "D");
		auto parameters = cart_prod(beta, J, D);
		
		write_logfile(registry, beta);
		RegularLatticeLoop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
	}




	else if (startswith(ham, "AshkinTeller"))
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto K    = registry.Get<std::vector<double> >("mc.ini", ham, "K");
		auto parameters = cart_prod(beta, J, K);
		
		write_logfile(registry, beta);
		RegularLatticeLoop<AshkinTeller<int>>(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "XXZAntiferro")
	{
		auto beta     = registry.Get<std::vector<double>>("mc.ini", ham, "beta");
		auto extfield = registry.Get<std::vector<double>>("mc.ini", ham, "extfield");
		auto aniso    = registry.Get<std::vector<double>>("mc.ini", ham, "aniso");
		auto parameters = cart_prod(beta, aniso, extfield);
		
		write_logfile(registry, beta);
		RegularLatticeLoop<XXZAntiferro<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "XXZAntiferroSingleAniso")
	{
		auto beta        = registry.Get<std::vector<double>>("mc.ini", ham, "beta");
		auto extfield    = registry.Get<std::vector<double>>("mc.ini", ham, "extfield");
		auto aniso       = registry.Get<std::vector<double>>("mc.ini", ham, "aniso");
		auto singleaniso = registry.Get<std::vector<double>>("mc.ini", ham, "singleaniso");
		auto parameters = cart_prod(beta, extfield, aniso, singleaniso);

		write_logfile(registry, extfield);
		RegularLatticeLoop<XXZAntiferroSingleAniso<double,double> >(registry, outbasedir, parameters, xxzfilter);
	}




	else if (startswith(ham, "Bimodal-Ising-EdwardsAnderson"))
	{
		const auto ham        = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc.ini", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc.ini", ham, "L" );
        int nthreads = 0;
        try {
            nthreads = registry.Get<int>("mc.ini", "General", "threads_per_node" );            
            }
        catch (const Registry_Key_not_found_Exception&) {
            std::cout<<"threads_per_node not set -> automatic"<<std::endl;
        }

		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);

        typedef decltype(finalize_parameter_triple(std::declval<std::tuple<int, int> >() ,std::declval<MARQOV::Config>(), hp)) PPType;
        typedef EdwardsAndersonIsing<int> Hamiltonian;
        typedef RegularRandomBond<BimodalPDF> Lattice;
        typename GetSchedulerType<Hamiltonian, Lattice, typename PPType::value_type>::MarqovScheduler sched(1, nthreads);

		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
	     	MARQOV::Config mp(outpath);
	     	mp.setnsweeps(15);
			mp.setncluster(0);
			mp.setwarmupsteps(300);
			mp.setgameloopsteps(300);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);
			// schedule simulations
			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter_triple);
//		 	Loop<Hamiltonian, RegularRandomBond<BimodalPDF> >(rparams, defaultfilter_triple);
		}
		// run!
		sched.start();
	}




	else if (startswith(ham, "Gaussian-Ising-EdwardsAnderson"))
	{
		const auto ham        = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc.ini", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc.ini", ham, "L" );
        int nthreads = 0;
        try {
            nthreads = registry.Get<int>("mc.ini", "General", "threads_per_node" );            
            }
        catch (const Registry_Key_not_found_Exception&) {
            std::cout<<"threads_per_node not set -> automatic"<<std::endl;
        }

		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);
        typedef decltype(finalize_parameter_triple(std::declval<std::tuple<int, int> >() ,std::declval<MARQOV::Config>(), hp)) PPType;
        typedef EdwardsAndersonIsing<int> Hamiltonian;
        typedef RegularRandomBond<GaussianPDF> Lattice;
        typename GetSchedulerType<Hamiltonian, Lattice, typename PPType::value_type>::MarqovScheduler sched(1, nthreads);
	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
        	MARQOV::Config mp(outpath);
        	mp.setnsweeps(50);
			mp.setncluster(0);
			mp.setwarmupsteps(100);
			mp.setgameloopsteps(1000);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

			// schedule simulations
			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter_triple);
//		 	Loop< EdwardsAndersonIsing<int>, RegularRandomBond<GaussianPDF>>(rparams, defaultfilter_triple);
		}
		sched.start(); // run!
	}




	else if (ham == "IsingCC")
	{
		const auto ham        = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc.ini", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc.ini", ham, "L" );
        int nthreads = 0;
        try {
            nthreads = registry.Get<int>("mc.ini", "General", "threads_per_node" );            
            }
        catch (const Registry_Key_not_found_Exception&) {
            std::cout<<"threads_per_node not set -> automatic"<<std::endl;
        }
		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);
        
        typedef decltype(finalize_parameter_triple(std::declval<std::tuple<int, int> >() ,std::declval<MARQOV::Config>(), hp)) PPType;
        typedef Ising<int> Hamiltonian;
        typedef ConstantCoordinationLattice<Poissonian> Lattice;
        typename GetSchedulerType<Hamiltonian, Lattice, typename PPType::value_type>::MarqovScheduler sched(1, nthreads);

		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
        	MARQOV::Config mp(outpath);
        	mp.setnsweeps(5);
			mp.setncluster(15);
			mp.setwarmupsteps(500);
			mp.setgameloopsteps(1500);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

            for (auto p: rparams)
                sched.createSimfromParameter(p, defaultfilter_triple);
			// perform simulations
		 	// Loop<Ising<int>, ConstantCoordinationLattice<Poissonian>>(rparams, defaultfilter_triple);
		}
		sched.start();
	}
	else if (ham == "IrregularIsing1")
	{
		// construct irregular lattice and pass it to MARQOV as a reference

		const int L = 32;
		const int dim = 2;
		ConstantCoordinationLattice<Poissonian> ccl(L, dim);

		// prepare output
		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		makeDir(outpath);

		// Hamiltonian parameters
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		std::vector<double> myj = {-1.0};
		auto hp = cart_prod(beta, myj);

		// Monte Carlo parameters
		MARQOV::Config mp(outpath);
		mp.setrepid(1);
		mp.setnsweeps(L);
		mp.setncluster(10);

		auto params = finalize_parameter_pair(mp, hp);
		// partially apply filter
		auto f = [&ccl](auto p){return defaultfilter(ccl, p);};

		typedef typename decltype(params)::value_type PPType;
		typedef Ising<int> Hamiltonian;
		typedef ConstantCoordinationLattice<Poissonian> Lattice;
		typename GetSchedulerType<Hamiltonian, Lattice, PPType>::MarqovScheduler sched(1);

		// schedule simulations
		for (auto p: params) sched.createSimfromParameter(p, f);
//		Loop<Ising<int>, ConstantCoordinationLattice<Poissonian>>(params, f);

		// run!
		sched.start();
	}

	else if (ham == "BlumeCapelBipartite")
	{
		// import parameters
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto DA   = registry.Get<std::vector<double> >("mc.ini", ham, "DA");
		auto DB   = registry.Get<std::vector<double> >("mc.ini", ham, "DB");
		auto parameters = cart_prod(beta, J, DA, DB);
        int nthreads = 0;
        try {
            nthreads = registry.Get<int>("mc.ini", "General", "threads_per_node" );            
            }
        catch (const Registry_Key_not_found_Exception&) {
            std::cout<<"threads_per_node not set -> automatic"<<std::endl;
        }

		const auto name      = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
		      auto nreplicas = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
		const auto nL 	      = registry.Get<std::vector<int>>("mc.ini", name, "L" );
		const auto dim 	 = registry.Get<int>("mc.ini", name, "dim" );

		write_logfile(registry, beta);
		
		std::vector<SimpleBipartite> latts;
        typedef decltype(finalize_parameter_pair(std::declval<MARQOV::Config>(), parameters)) PPType;
        typename GetSchedulerType<BlumeCapelBipartite<int>, SimpleBipartite, typename PPType::value_type>::MarqovScheduler sched(1, nthreads);
		// set up replicas
		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
		
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare. Extend lifetime of lattices.
			int L = nL[j];
			latts.emplace_back(L, dim);
		}
		
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare
			int L = nL[j];

			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			MARQOV::Config mp(outpath);
			makeDir(mp.outpath);

			mp.setnsweeps(2);
			mp.setncluster(int(L/2));
			mp.setwarmupsteps(500);
			mp.setgameloopsteps(2500);

			// set up parameter space
			auto params = finalize_parameter_pair(mp, parameters);
			auto rparams = replicator_pair(params, nreplicas[j]);

			// lattice
//			SimpleBipartite latt(L, dim);
			SimpleBipartite& latt = latts[j]; 

			// partially apply filter
	 		auto f = [&latt, &outbasedir, L](auto p){return defaultfilter(latt, p);};

			// schedule
			for(auto p : rparams) sched.createSimfromParameter(p, f);
//	 		Loop<BlumeCapelBipartite<int>, SimpleBipartite>(rparams, f);
		}
		sched.start(); // run!
	}
}



int main()
{
    std::cout<<"MARQOV Copyright (C) 2020-2021, The MARQOV Project contributors"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;

	// read config files
	RegistryDB registry("../src/config", "ini");

	// remove old output and prepare new one
	std::string outbasedir = registry.Get<std::string>("mc.ini", "IO", "outdir" );
	std::string logbasedir = registry.Get<std::string>("mc.ini", "IO", "logdir" );

	//FIXME: NEVER DELETE USER DATA
	std::string command;
	command = "rm -r " + outbasedir;
	system(command.c_str());
	command = "rm -r " + logbasedir;
	system(command.c_str());

	makeDir(outbasedir);
	makeDir(logbasedir);

	selectsim(registry, outbasedir, logbasedir);
}
