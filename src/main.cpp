#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <iomanip>

#define SSH_2D  // dimension switch

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

#include "rndwrapper.h"
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
#include "hamiltonian/Ssh.h"
#include "hamiltonian/BlumeCapelBipartite.h"

using namespace MARQOV;


// ---------------------------------------

std::string selectsim_startup(RegistryDB& registry)
{
	const auto ham        = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
	const auto dim 	  = registry.Get<int>("mc.ini", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>("mc.ini", ham, "rep" );
	const auto nreplicass = registry.Get<std::string>("mc.ini", ham, "rep" );
	const auto nL  	  = registry.Get<std::vector<int>>("mc.ini", ham, "L" );
	const auto nLs 	  = registry.Get<std::string>("mc.ini", ham, "L" );

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

     auto beta   = registry.Get<std::vector<double> >("mc.ini", ham, "betaMC");
     auto dtau   = registry.Get<std::vector<double> >("mc.ini", ham, "dtau");
     auto m      = registry.Get<std::vector<double> >("mc.ini", ham, "m");
     auto k      = registry.Get<std::vector<double> >("mc.ini", ham, "k");
     auto g      = registry.Get<std::vector<double> >("mc.ini", ham, "g");
     auto mu     = registry.Get<std::vector<double> >("mc.ini", ham, "mu");

     const auto name      = registry.Get<std::string>("mc.ini", "General", "Hamiltonian" );
           auto nreplicas = registry.Get<std::vector<int>>("mc.ini", name, "rep" );
     const auto nL        = registry.Get<std::vector<int>>("mc.ini", name, "L" );
     const auto nLtime    = registry.Get<std::vector<int>>("mc.ini", name, "Ltime" );
     const auto dim       = registry.Get<int>("mc.ini", name, "dim" );


	// we need "L" and "Ltime" as explicit parameters in the Hamiltonian
	// which requires some gymnastics ...
	std::vector<int> dummy = {0};
     auto hp = cart_prod(beta, m, k, g, mu, dtau, dummy, dummy);


	// prepare lattices
	std::vector<std::vector<SSHLattice>> latts;
	for (std::size_t j=0; j<nL.size(); j++)
	{
		std::vector<SSHLattice> lat;

     	for (std::size_t jj=0; jj<nLtime.size(); jj++)
		{
			const int L     = nL[j];
			const int Ltime = nLtime[jj];


			lat.emplace_back(L, Ltime, dim);
		}
		latts.emplace_back(lat);
	}

	typedef decltype(finalize_parameter_pair(std::declval<MARQOV::Config>(), hp)) PPType; 
	typename GetSchedulerType<SSH<double>, SSHLattice, typename PPType::value_type>::MarqovScheduler sched(1);

     for (std::size_t j=0; j<nL.size(); j++)
	{
     	for (std::size_t jj=0; jj<nLtime.size(); jj++)
		{

			const int L     = nL[j];
			const int Ltime = nLtime[jj];
     	     cout << endl << "L_space = " << L << "\t" << "L_time = " << Ltime << endl << endl;


			for (std::size_t i=0; i<hp.size(); i++) std::get<6>(hp[i]) = Ltime;
			for (std::size_t i=0; i<hp.size(); i++) std::get<7>(hp[i]) = L;

     	     std::string outpath = outbasedir+"/"+std::to_string(L)+"/";

     	     MARQOV::Config mp(outpath);
     	     mp.setnsweeps(10);
     	     mp.setncluster(0);
     	     mp.setwarmupsteps(0);
     	     mp.setgameloopsteps(300);

     	     makeDir(mp.outpath);


     	     // set up parameters
     	     auto rparams = finalize_parameter_pair(mp, hp);
			auto params = replicator_pair(rparams, nreplicas[0]);

     	     SSHLattice& latt = latts[j][jj];
     	     auto f = [&latt, &outbasedir, L](auto p){return sshfilter(latt, p);}; //partially apply filter
			for(auto p : params)
				sched.createSimfromParameter(p, f);
     	}
     }
	sched.start();

}







int main()
{

	// read config files
	RegistryDB registry("../src/config", "ini");

	// remove old output and prepare new one
	std::string outbasedir = registry.Get<std::string>("mc.ini", "IO", "outdir" );
	std::string logbasedir = registry.Get<std::string>("mc.ini", "IO", "logdir" );

	//FIXME: NEVER DELETE USER DATA
	std::string command;
	command = "rm -r " + outbasedir;
	system(command.c_str());
//	command = "rm -r " + logbasedir;
//	system(command.c_str());

	makeDir(outbasedir);
//	makeDir(logbasedir);
	
	selectsim(registry, outbasedir, logbasedir);
}
