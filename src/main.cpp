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
#include "embedder.h"
#include "marqovscheduler.h"
#include "util.h"

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
#include "hamiltonian/AshkinTeller.h"

using namespace MARQOV;



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



	if (ham == "Ising")
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto parameters = cart_prod(beta, J);

 		RegularLatticeLoop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
	}
//
//
//	else if (ham == "AshkinTeller")
//	{
//		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
//		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
//		auto K    = registry.Get<std::vector<double> >("mc.ini", ham, "K");
//		auto parameters = cart_prod(beta, J, K);
//
//		RegularLatticeLoop<AshkinTeller<int> >(registry, outbasedir, parameters, defaultfilter);
//	}
//
//
	else if (ham == "Heisenberg")
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto parameters = cart_prod(beta, J);

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
		for (std::size_t i=0; i<parameters.size(); i++) 
			std::get<1>(parameters[i]) = std::get<0>(parameters[i]);
		
		RegularLatticeLoop<Phi4<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}




	else if (ham == "BlumeCapel")
	{
		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
		auto D    = registry.Get<std::vector<double> >("mc.ini", ham, "D");
		auto parameters = cart_prod(beta, J, D);
		
		RegularLatticeLoop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
	}




//	else if (startswith(ham, "AshkinTeller"))
//	{
//		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
//		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
//		auto K    = registry.Get<std::vector<double> >("mc.ini", ham, "K");
//		auto parameters = cart_prod(beta, J, K);
//		
//		RegularLatticeLoop<AshkinTeller<int>>(registry, outbasedir, parameters, defaultfilter);
//	}




	else if (ham == "XXZAntiferro")
	{
		auto beta     = registry.Get<std::vector<double>>("mc.ini", ham, "beta");
		auto extfield = registry.Get<std::vector<double>>("mc.ini", ham, "extfield");
		auto aniso    = registry.Get<std::vector<double>>("mc.ini", ham, "aniso");
		auto parameters = cart_prod(beta, aniso, extfield);
		
		RegularLatticeLoop<XXZAntiferro<double>>(registry, outbasedir, parameters, defaultfilter);
	}




//	else if (ham == "XXZAntiferroSingleAniso")
//	{
//		auto beta        = registry.Get<std::vector<double>>("mc.ini", ham, "beta");
//		auto extfield    = registry.Get<std::vector<double>>("mc.ini", ham, "extfield");
//		auto aniso       = registry.Get<std::vector<double>>("mc.ini", ham, "aniso");
//		auto singleaniso = registry.Get<std::vector<double>>("mc.ini", ham, "singleaniso");
//		auto parameters = cart_prod(beta, extfield, aniso, singleaniso);
//
//		RegularLatticeLoop<XXZAntiferroSingleAniso<double,double> >(registry, outbasedir, parameters, xxzfilter);
//	}
//
//
//
//	else if (startswith(ham, "EdwardsAnderson-Ising"))
//	{
//		// Parameters
//		const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
//		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
//		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
//		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
//	
//	
//		// Number of threads
//		int nthreads = 0;
//		try 
//		{
//			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
//		}
//		catch (const Registry_Key_not_found_Exception&) 
//		{
//			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
//		}
//
//		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
//
//		// Physical parameters
//		auto beta = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "beta");
//		auto J    = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "J");
//		auto hp = cart_prod(beta, J);
//        
//
//		// Typedefs
//		typedef EdwardsAndersonIsing<int> Hamiltonian;
//
//		typedef RegularRandomBond<GaussianPDF> Lattice;
//        	//typedef RegularRandomBond<BimodalPDF> Lattice;
//
//		typedef decltype(finalize_parameter_triple(std::declval<std::tuple<int, int> >() ,std::declval<MARQOV::Config>(), hp)) ParameterTripleType;
//		typedef typename ParameterTripleType::value_type ParameterType;
//		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;
//
//		SchedulerType sched(1);
//
//		// Lattice size loop
//		for (std::size_t j=0; j<nL.size(); j++)
//		{
//			// prepare output
//			int L = nL[j];
//			cout << endl << "L = " << L << endl << endl;
//			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
//			makeDir(outpath);
//	
//			// Monte Carlo parameters
//			MARQOV::Config mp(outpath);
//			mp.setnsweeps(50);
//			mp.setncluster(0);
//			mp.setwarmupsteps(200);
//			mp.setgameloopsteps(1000);
//
//			// lattice parameters
//			auto lp = std::make_tuple(L,dim);
//
//			// form parameter triple and replicate
//			auto params  = finalize_parameter_triple(lp, mp, hp);
//			auto rparams = replicator(params, nreplicas[j]);
//
//			// schedule simulations
//			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter_triple);
//		}
//		sched.start(); // run!
//	}
//
//
//
//
//	else if (ham == "IsingCC")
//	{
//		// Parameters
//		const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
//		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
//		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
//		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
//	
//	
//		// Number of threads
//		int nthreads = 0;
//		try 
//		{
//			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
//		}
//		catch (const Registry_Key_not_found_Exception&) 
//		{
//			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
//		}
//
//
//		// Replicas
//		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
//
//		// Physical parameters
//		auto beta = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "beta");
//		auto J    = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "J");
//		auto hp = cart_prod(beta, J);
//        
//
//		// Typedefs
//		typedef Ising<int> Hamiltonian;
//		typedef ConstantCoordinationLattice<Poissonian> Lattice;
//		typedef decltype(finalize_parameter_triple(std::declval<std::tuple<int, int> >() ,std::declval<MARQOV::Config>(), hp)) ParameterTripleType;
//		typedef typename ParameterTripleType::value_type ParameterType;
//		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;
//
//
//		// Lattice size loop
//		for (std::size_t j=0; j<nL.size(); j++)
//		{
//			// init scheduler
//			SchedulerType sched(1, nthreads);
//
//			// prepare output
//			int L = nL[j];
//			cout << endl << "L = " << L << endl << endl;
//			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
//			makeDir(outpath);
//	
//			// Monte Carlo parameters
//			MARQOV::Config mp(outpath);
//			mp.setnsweeps(5);
//			mp.setncluster(15);
//			mp.setwarmupsteps(500);
//			mp.setgameloopsteps(1500);
//
//			// lattice parameters
//			auto lp = std::make_tuple(L,dim);
//
//			// form parameter triple and replicate
//			auto params  = finalize_parameter_triple(lp, mp, hp);
//			auto rparams = replicator(params, nreplicas[j]);
//
//			// feed scheduler
//			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter_triple);
//
//			// run!
//			sched.start();
//		}
//	}
//
//
//
//
//
//
//
//
//	else if (ham == "BlumeCapelBipartite")
//	{
//		// Parameters
//		const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
//		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
//		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
//		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
//	
//	
//		// Number of threads
//		int nthreads = 0;
//		try 
//		{
//			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
//		}
//		catch (const Registry_Key_not_found_Exception&) 
//		{
//			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
//		}
//
//
//		// Replicas
//		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
//
//
//		// import parameters
//		auto beta = registry.Get<std::vector<double> >("mc.ini", ham, "beta");
//		auto J    = registry.Get<std::vector<double> >("mc.ini", ham, "J");
//		auto DA   = registry.Get<std::vector<double> >("mc.ini", ham, "DA");
//		auto DB   = registry.Get<std::vector<double> >("mc.ini", ham, "DB");
//		auto hp = cart_prod(beta, J, DA, DB);
//
//		typedef BlumeCapelBipartite<int> Hamiltonian;
//		typedef SimpleBipartite Lattice;
//
//
//		typedef decltype(finalize_parameter_pair(std::declval<MARQOV::Config>(), hp)) ParameterPairType;
//		typedef typename ParameterPairType::value_type ParameterType;
//		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType>::MarqovScheduler SchedulerType;
//		
//		
//		// Prepare Geometry
//		std::vector<SimpleBipartite> latts;
//		for (std::size_t j=0; j<nL.size(); j++) latts.emplace_back(nL[j], dim);
//	    
//	
//		// Init Scheduler
//		SchedulerType sched(1, nthreads);
//	    
//	
//		// Lattice size loop
//		for (std::size_t j=0; j<nL.size(); j++)
//		{
//			// prepare
//			int L = nL[j];
//			cout << endl << "L = " << L << endl << endl;
//	
//			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
//	
//			MARQOV::Config mp(outpath);
//			mp.setnsweeps(10);
//			mp.setncluster(int(L/2));
//			mp.setwarmupsteps(200);
//			mp.setgameloopsteps(500);
//			
//			makeDir(mp.outpath);
//			
//			auto params = finalize_parameter_pair(mp, hp);
//			auto rparams = replicator_pair(params, nreplicas[j]);
//			
//			// set up and execute        
//			Lattice& latt = latts[j];
//			auto f = [&latt](auto p){return defaultfilter(latt, p);}; //partially apply filter
//	
//			// feed the scheduler
//			for(auto p: rparams) sched.createSimfromParameter(p, f);
//		}
//		sched.start();
//	}
}






int main()
{



//	Heisenberg<double, double> ham = Heisenberg<double, double>(1);
//	std::array<double, 3> arr = {1,2,3};
//	std::vector<int> nbrs = {11,111,1111};
//	auto k = wolff_embedding<Heisenberg<double, double>, std::array<double, 3>, std::vector<int>>(ham, arr, nbrs);
//	cout << k << endl;
//	exit(0);


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
