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
#include "hamiltonian/BlumeCapelBipartite.h"
#include "hamiltonian/AshkinTeller.h"

using namespace MARQOV;

/** Find out if a string starts with something.
 * @param longword we search in this string
 * @param shortword we look for this
 * @return truen if longword strarts with shortword, else false
 */
bool startswith(const std::string& longword, const std::string& shortword) noexcept
{
    return longword.find(shortword) == 0;
}

void checkreplicaconfig(int nr, int nL)
{
    if ((nr != nL) && (nr != 1)) throw std::invalid_argument("[MARQOV] Invalid replica configuration!");
}

void tidyupoldsims(const std::string& outbasedir)
{
	// delete previous output // fixme: don't do that by default!
#ifdef MPIMARQOV
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if(myrank == 0)
    {
#endif
        std::cout<<"[MARQOV::main] Erasing previous data!!!!!!!!!!!!!!!"<<std::endl;
        std::string command = "rm -r " + outbasedir;
        system(command.c_str());
        makeDir(outbasedir);
#ifdef MPIMARQOV
    }
#endif
}

void printInfoandcheckreplicaconfig(RegistryDB& registry, const std::string& ham)
{
    const auto dim 	      = registry.Get<int>(ham+".ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>(ham+".ini", ham, "rep" );
	const auto nreplicas_str = registry.Get<std::string>(ham+".ini", ham, "rep" );
	const auto nL  	      = registry.Get<std::vector<int>>(ham+".ini", ham, "L" );
	const auto nLs 	      = registry.Get<std::string>(ham+".ini", ham, "L" );

	cout << endl;
	cout << "Hamiltonian: \t" << ham << endl;
	cout << "Dimension: \t" << dim << endl;
	cout << "Lattice sizes:\t" << nLs << endl;
	cout << "Replicas:\t" << nreplicas_str << endl;
    
    checkreplicaconfig(nreplicas.size(), nL.size());
}

void createcfgfooter(std::ostream& os, int nmetro, double nclusteramp, int nclusterexp, int warmupsteps, int measuresteps)
{
    os<<"[MC]\n"<<"nmetro = "<<nmetro<<"\nnclusteramp = "<<nclusteramp<<"\nnclusterexp = "<<nclusterexp<<"\nwarmupsteps = "<<warmupsteps<<"\nmeasuresteps = "<<measuresteps<<"\n";
    os<<"[IO]\n"<<"outdir = ../out\n"<<"[END]"<<std::endl;
}

void scheduleIsing(RegistryDB& registry)
{
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>("Ising.ini", "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find Ising config! Generating new one in ./config/Ising.ini"<<std::endl;
        ofstream ising("./config/Ising.ini");
        ising<<"[Ising]\n"<<"L = 10\n"<<"rep = 2\n"<<"dim = 2\n"<<"beta = 0.3, 0.4\n"<<"J = -1.0\n\n";
        createcfgfooter(ising, 10, 25, 0, 30, 30);
        registry.init("./config");
        outbasedir = registry.Get<std::string>("Ising.ini", "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "Ising");
    auto beta = registry.Get<std::vector<double> >("Ising.ini", "Ising", "beta");
    auto J    = registry.Get<std::vector<double> >("Ising.ini", "Ising", "J");
    auto parameters = cart_prod(beta, J);
    RegularLatticeLoop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
}

void schedulePotts(RegistryDB& registry)
{
    std::string ham = registry.Get<std::string>("select.ini", "General", "Hamiltonian" );
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(ham+".ini", "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        int q = std::stoi(ham.substr(5));
        std::cout<<"[MARQOV] Unable to find"<<ham<<" config! Generating new one in ./config/"<<ham<<".ini"<<std::endl;
        ofstream cfg("./config/" + ham + ".ini");
        cfg<<"["<<ham<<"]\n"<<"q = "<<q<<'\n'<<"L = 10\n"<<"rep = 5\n"<<"dim = 2\n"<<"beta = 1.05, 1.07, 1.09\n"<<"J = -1.0\n\n";
        createcfgfooter(cfg, 2, 10, 0, 200, 200);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(ham + ".ini", "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, ham);
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

void scheduleAshkinTeller(RegistryDB& registry)
{
    const std::string fn{"Ashkin-Teller.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find Ashkin-Teller config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[AshkinTeller]\n"<<"L = 12\n"<<"rep = 10\n"<<"dim = 2\n"<<"beta = 0.312, 0.313\n"<<"J = -1.0\n"<<"K = 0.3\n\n";
        createcfgfooter(cfg, 2, 25, 0, 500, 10000);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "AshkinTeller");
    auto beta = registry.Get<std::vector<double> >(fn, "AshkinTeller", "beta");
    auto J    = registry.Get<std::vector<double> >(fn, "AshkinTeller", "J");
    auto K    = registry.Get<std::vector<double> >(fn, "AshkinTeller", "K");
    auto parameters = cart_prod(beta, J, K);

    RegularLatticeLoop<AshkinTeller>(registry, outbasedir, parameters, defaultfilter);
}

void scheduleHeisenberg(RegistryDB& registry)
{
    const std::string fn{"Heisenberg.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find Heisenberg config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[Heisenberg]\n"<<"L = 8,10,12\n"<<"rep = 4\n"<<"dim = 3\n"<<"beta = 0.67,0.68\n"<<"J = -1.0\n\n";
        createcfgfooter(cfg, 2, 0.5, 1, 500, 5000);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "Heisenberg");
    auto beta = registry.Get<std::vector<double> >(fn, "Heisenberg", "beta");
    auto J    = registry.Get<std::vector<double> >(fn, "Heisenberg", "J");
    auto parameters = cart_prod(beta, J);
    
    RegularLatticeLoop<Heisenberg<double, double> >(registry, outbasedir, parameters, defaultfilter);
}

void schedulePhi4(RegistryDB& registry)
{
    const std::string fn{"Phi4.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find Phi4 config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[Phi4]\n"<<"L = 6,8,12\n"<<"rep = 10\n"<<"dim = 3\n"<<"beta = 0.681, 0.684, 0.687\n"<<"lambda = 4.5\n"<<"mass = 1.0\n\n";
        createcfgfooter(cfg, 5, 0.5, 1, 20, 100);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "Phi4");
    auto beta   = registry.Get<std::vector<double> >(fn, "Phi4", "beta");
    auto lambda = registry.Get<std::vector<double> >(fn, "Phi4", "lambda");
    auto mass   = registry.Get<std::vector<double> >(fn, "Phi4", "mass");
		
		// we need "beta" as an explicit parameter in the Hamiltonian
		// this requires some gymnastics ...
		std::vector<double> dummy = {0.0};
		auto parameters = cart_prod(beta, dummy, lambda, mass);
		for (std::size_t i=0; i<parameters.size(); i++) 
			std::get<1>(parameters[i]) = std::get<0>(parameters[i]);
		
		RegularLatticeLoop<Phi4<double, double> >(registry, outbasedir, parameters, defaultfilter);
}

void scheduleBlumeCapel(RegistryDB& registry)
{
    const std::string fn{"BlumeCapel.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find BlumeCapel config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[BlumeCapel]\n"<<"L = 8,12\n"<<"rep = 2\n"<<"dim = 3\n"<<"beta = 0.32,0.33\n"<<"J = -1\n"<<"D = 0.655\n\n";
        createcfgfooter(cfg, 3, 0.5, 1, 100, 1000);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "BlumeCapel");
    auto beta = registry.Get<std::vector<double> >(fn, "BlumeCapel", "beta");
    auto J    = registry.Get<std::vector<double> >(fn, "BlumeCapel", "J");
    auto D    = registry.Get<std::vector<double> >(fn, "BlumeCapel", "D");
    auto parameters = cart_prod(beta, J, D);
    
    RegularLatticeLoop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
}

void scheduleBlumeEmeryGriffiths(RegistryDB& registry)
{
    const std::string fn{"BlumeEmeryGriffiths.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find BlumeEmeryGriffiths config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[BlumeEmeryGriffiths]\n"<<"L = 8,12\n"<<"rep = 4\n"<<"dim = 2\n"<<"beta = 0.42\n"<<"J = -1\n"<<"D = 1\n"<<"K = 1\n\n";
        createcfgfooter(cfg, 2, 0.5, 1, 300, 3000);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "BlumeEmeryGriffiths");
		auto beta = registry.Get<std::vector<double> >(fn, "BlumeEmeryGriffiths", "beta");
		auto J    = registry.Get<std::vector<double> >(fn, "BlumeEmeryGriffiths", "J");
		auto D    = registry.Get<std::vector<double> >(fn, "BlumeEmeryGriffiths", "D");
		auto K    = registry.Get<std::vector<double> >(fn, "BlumeEmeryGriffiths", "K");
		auto parameters = cart_prod(beta, J, D, K);

		RegularLatticeLoop<BlumeEmeryGriffiths<int>>(registry, outbasedir, parameters, defaultfilter);
}

void scheduleXXZAntiferro(RegistryDB& registry)
{
    const std::string fn{"XXZAntiferro.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find XXZAntiferro config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[XXZAntiferro]\n"<<"L = 8,10,12\n"<<"rep = 10\n"<<"dim = 3\n"<<"beta = 0.92,0.93\n"<<"J = -1\n"<<"extfield = 1\n"<<"aniso = 0.8\n\n";
        createcfgfooter(cfg, 3, 0.0, 1, 500, 1500);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "XXZAntiferro");
		auto beta     = registry.Get<std::vector<double>>("XXZAntiferro.ini", "XXZAntiferro", "beta");
		auto extfield = registry.Get<std::vector<double>>("XXZAntiferro.ini", "XXZAntiferro", "extfield");
		auto aniso    = registry.Get<std::vector<double>>("XXZAntiferro.ini", "XXZAntiferro", "aniso");
		auto parameters = cart_prod(beta, aniso, extfield);

		RegularLatticeLoop<XXZAntiferro<double>>(registry, outbasedir, parameters, defaultfilter);
}

void scheduleXXZAntiferroSingleAniso(RegistryDB& registry)
{
    const std::string fn{"XXZAntiferroSingleAniso.ini"};
    std::string outbasedir;
    try
    {
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    catch(Registry_cfgfile_not_found_Exception& e) 
    {
        std::cout<<"[MARQOV] Unable to find XXZAntiferroSingleAniso config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg<<"[XXZAntiferroSingleAniso]\n"<<"L = 8,10,12\n"<<"rep = 10\n"<<"dim = 3\n"<<"beta = 0.92,0.93\n"<<"J = -1\n"<<"extfield = 1\n"<<"aniso = 0.8\n"<<"singleaniso = 0.4\n\n";
        createcfgfooter(cfg, 3, 0.0, 1, 500, 1500);
        registry.init("./config");
        outbasedir = registry.Get<std::string>(fn, "IO", "outdir" );
    }
    tidyupoldsims(outbasedir);
    printInfoandcheckreplicaconfig(registry, "XXZAntiferroSingleAniso");
        auto beta        = registry.Get<std::vector<double>>(fn, "XXZAntiferroSingleAniso", "beta");
		auto extfield    = registry.Get<std::vector<double>>(fn, "XXZAntiferroSingleAniso", "extfield");
		auto aniso       = registry.Get<std::vector<double>>(fn, "XXZAntiferroSingleAniso", "aniso");
		auto singleaniso = registry.Get<std::vector<double>>(fn, "XXZAntiferroSingleAniso", "singleaniso");
		auto parameters = cart_prod(beta, extfield, aniso, singleaniso);

		RegularLatticeLoop<XXZAntiferroSingleAniso<double>>(registry, outbasedir, parameters, xxzfilter);
}

void scheduleEdwardsAndersonIsing(RegistryDB& registry)
{
// Parameters
    const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
    std::string outbasedir = registry.Get<std::string>(name + ".ini", "IO", "outdir" );
    tidyupoldsims(outbasedir);
		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
        printInfoandcheckreplicaconfig(registry, name);	
	
		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}

		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		// Physical parameters
		auto beta = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "J");
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

			// form parameter triple and replicate
			auto params  = finalize_parameter(std::make_tuple(L, dim)/*lattice parameters*/, mp, hp);//this particular form is required to happify PGI-19.10
            auto rparams = replicator(params, nreplicas[j]);

			// schedule simulations
 			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);
		}
 		sched.start(); // run!
}

void scheduleIsingCC(RegistryDB& registry)
{
		// Parameters
    const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
    std::string outbasedir = registry.Get<std::string>(name + ".ini", "IO", "outdir" );
    tidyupoldsims(outbasedir);
		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
        printInfoandcheckreplicaconfig(registry, name);

		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}


		// Replicas
		if (nreplicas.size() == 1) { for (decltype(nL.size()) i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		// Physical parameters
		auto beta = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", "IsingCC", "J");
		auto hp = cart_prod(beta, J);

		// Typedefs
		typedef Ising<int> Hamiltonian;
		typedef ConstantCoordinationLattice<Poissonian> Lattice;

        typedef std::tuple<std::tuple<int, int>, MARQOV::Config, typename decltype(hp)::value_type > ParameterType;
		typedef typename GetSchedulerType<Hamiltonian, Lattice, ParameterType, std::knuth_b>::MarqovScheduler SchedulerType;


		// Lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// init scheduler
			SchedulerType sched(1, nthreads);

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

			// feed scheduler
			for (auto p: rparams) sched.createSimfromParameter(p, defaultfilter);

			// run!
			sched.start();
		}
}

void scheduleBlumeCapelBiPartite(RegistryDB& registry)
{
    // Parameters
    const auto name = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
    std::string outbasedir = registry.Get<std::string>(name + ".ini", "IO", "outdir" );
    tidyupoldsims(outbasedir);
		auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
		const auto nL   = registry.Get<std::vector<int>>("mc.ini", name, "L" );
		const auto dim  = registry.Get<int>("mc.ini", name, "dim" );
        printInfoandcheckreplicaconfig(registry, name);

		// Number of threads
		int nthreads = 0;
		try 
		{
			nthreads = registry.template Get<int>("mc.ini", "General", "threads_per_node" );
		}
		catch (const Registry_Key_not_found_Exception&) 
		{
			std::cout<<"threads_per_node not set -> automatic"<<std::endl;
		}


		// Replicas
		if (nreplicas.size() == 1) { for (std::size_t i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }


		// import parameters
		auto beta = registry.Get<std::vector<double> >("mc.ini", "BlumeCapelBipartite", "beta");
		auto J    = registry.Get<std::vector<double> >("mc.ini", "BlumeCapelBipartite", "J");
		auto DA   = registry.Get<std::vector<double> >("mc.ini", "BlumeCapelBipartite", "DA");
		auto DB   = registry.Get<std::vector<double> >("mc.ini", "BlumeCapelBipartite", "DB");
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

/** Select the respective simulation.
 * 
 * @param registry The registry object that we will use
 */
void selectsim(RegistryDB& registry)
{
	const auto ham = registry.Get<std::string>("select.ini", "General", "Hamiltonian" );

	// ----------------- select simulation ------------------

	if (startswith(ham, "Ising"))
        scheduleIsing(registry);

	else if (startswith(ham,"Potts"))
        schedulePotts(registry);

	else if (ham == "AshkinTeller")
        scheduleAshkinTeller(registry);

	else if (ham == "Heisenberg")
        scheduleHeisenberg(registry);

	else if (ham == "Phi4")
        schedulePhi4(registry);

	else if (ham == "BlumeCapel")
        scheduleBlumeCapel(registry);

	else if (ham == "BlumeEmeryGriffiths")
        scheduleBlumeEmeryGriffiths(registry);

	else if (ham == "XXZAntiferro")
        scheduleXXZAntiferro(registry);

	else if (ham == "XXZAntiferroSingleAniso")
        scheduleXXZAntiferroSingleAniso(registry);

	else if (startswith(ham, "EdwardsAnderson-Ising"))
        scheduleEdwardsAndersonIsing(registry);

	else if (ham == "IsingCC")
        scheduleIsingCC(registry);

	else if (ham == "BlumeCapelBipartite")
        scheduleBlumeCapelBiPartite(registry);
}

int main(int argc, char* argv[])
{
    std::cout<<"MARQOV Copyright (C) 2020-2021, The MARQOV Project contributors"<<std::endl;
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<std::endl;
    std::cout<<"This is free software, and you are welcome to redistribute it under certain conditions."<<std::endl;
#ifdef MPIMARQOV
    int threadingsupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &threadingsupport);//FIXME: maybe we get by with one level less.
    if(threadingsupport < MPI_THREAD_SERIALIZED)
    {
        std::cout<<"[MARQOV::main] Couldn't initialize MPI! Requested threading level not supported."<<std::endl;
        return -1;
    }
#endif

	// read config files
	RegistryDB registry;
    try
    {
        registry.init("../src/config", "ini");
    }
    catch(Registry_Exception& re)
    {
        std::cout<<"[MARQOV::main] Configuration directory not found! Assuming you're starting the MARQOV demonstration binary for the first time!"<<std::endl;
        std::cout<<"WELCOME TO MARQOV!"<<std::endl;
        std::cout<<"[MARQOV::main] To get you going we will generate and populate a configuration directory locally under ./config"<<std::endl;
        makeDir("./config");
        const auto filename = std::string{"./config/select.ini"};
        if(!fileexists(filename))
        {
            std::ofstream select(filename);
            select<<"[General]"<<'\n'<<"Hamiltonian = Ising"<<'\n'<<"[END]"<<std::endl;
            registry.init("./config", "ini");
        }
        else
        {
            std::cout<<"[MARQOV::main] "<<filename<<" already exists, but is not usable. I would overwrite its content, hence I'm terminating now"<<std::endl;
            throw;
        }
    }
    int myrank = 0;
#ifdef MPIMARQOV
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
#endif
        selectsim(registry);
#ifdef MPIMARQOV
    }
    MPI_Finalize();
#endif
}
