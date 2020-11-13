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

#include "rndwrapper.h"
#include "helpers.h"
#include "geom/regular_lattice.h"
#include "geom/ssh_lattice.h"
#include "geom/grid.h"
#include "geom/neighbourclass.h"
#include "geom/io.h"
#include "vectorhelpers.h"
#include "cartprod.h"
#include "registry.h"
#include "systemtools.h"
#include "replicate.h"
#include "marqov.h"
#include "svmath.h"
#include "filters.h"

#include "marqovscheduler.h"
// Hamiltonians
#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"
#include "XXZAntiferroSingleAniso.h"
#include "AshkinTeller.h"
#include "EdwardsAndersonIsing.h"
#include "Ssh.h"
#include "BlumeCapelBipartite.h"

using namespace MARQOV;



//c++17 make_from_tuple from cppreference adapted for emplace
template <class Cont, class Tuple1, class Tuple2, std::size_t... I>
constexpr auto emplace_from_tuple_impl(Cont&& cont, Tuple1&& t1, MARQOV::MARQOVConfig&& mc, Tuple2&& t2, std::index_sequence<I...> )
{
	return cont.emplace_back(
		std::forward<Tuple1>(t1), 
		std::forward<MARQOV::MARQOVConfig>(mc), 
		std::get<I>(std::forward<Tuple2>(t2))...);
}
 
template <class Cont, class Tuple1, class Tuple2>
constexpr auto emplace_from_tuple(Cont&& cont, Tuple1&& t1, MARQOV::MARQOVConfig&& mc, Tuple2&& t2 )
{
	return emplace_from_tuple_impl(
		cont, 
		std::forward<Tuple1>(t1), 
		std::forward<MARQOV::MARQOVConfig>(mc), 
		std::forward<Tuple2>(t2),
		std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple2>>::value>{});
}

template <class Args1, class Args2, class T, class Callable>
void fillsims(const std::vector<Triple<Args1, MARQOV::MARQOVConfig, Args2> >& args, std::vector<T>& sims, Callable c)
{
	for(auto p : args)
	{
		auto t1 = c(p);
		emplace_from_tuple(sims, t1.first, std::forward<MARQOV::MARQOVConfig>(t1.second), t1.third);
	}
}

template <class Args, class T, class Callable>
void fillsims(const std::vector<std::pair<MARQOV::MARQOVConfig, Args>>& args, std::vector<T>& sims, Callable c)
{
	for(auto p : args)
	{
		auto t1 = c(p);
		emplace_from_tuple(sims, 
			std::forward<decltype(std::get<0>(t1))>(std::get<0>(t1)), 
			std::forward<MARQOV::MARQOVConfig>(std::get<1>(t1)), std::get<2>(t1));
	}
}

// ---------------------------------------


template<class ... Ts> struct sims_helper {};

template <class H,  class L, class HArgstuple, size_t... S>
struct sims_helper<H, L, HArgstuple, std::index_sequence<S...> >
{
	typedef decltype(makeMarqov<H>(std::declval<L>(),
	                       		 std::declval<MARQOVConfig>(),
	                       		 std::declval<typename std::tuple_element<S, HArgstuple>::type>()...
							 )) MarqovType;
};

template <class ... Ts>
struct sims_helper2 {};

template <class Hamiltonian, class Lattice, class LArgs, class HArgs>
struct sims_helper2<Hamiltonian, Lattice, Triple<LArgs, MARQOVConfig, HArgs> >
{
    typedef decltype(makeMarqov<Hamiltonian, Lattice>(std::declval<MARQOVConfig>(),  
    							  std::declval<std::pair<LArgs, HArgs>& >()
							  )) MarqovType;
};

template <class Hamiltonian, class Lattice, class HArgs>
struct sims_helper2<Hamiltonian, Lattice, std::pair<MARQOVConfig, HArgs> >
{
    static constexpr std::size_t tsize = std::tuple_size<typename std::remove_reference<HArgs>::type>::value;
    typedef std::make_index_sequence<tsize> HArgSequence;
    typedef typename sims_helper<Hamiltonian, Lattice, HArgs, HArgSequence>::MarqovType MarqovType;
};

// ---------------------------------------

template <class Hamiltonian, class Lattice, class Parameters, class Callable>
void Loop(const std::vector<Parameters>& params, Callable filter)
{
    typedef typename sims_helper2<Hamiltonian, Lattice, Parameters >::MarqovType MarqovType;  

    //create simulations
    std::vector<MarqovType> sims;
    sims.reserve(params.size());
    fillsims(params, sims, filter);
    
    Scheduler<typename decltype(sims)::value_type> sched(1);
    
    for (int i = 0; i < sims.size(); ++i)
     {
         sched.enqueuesim(sims[i]);
     }
    sched.start();
}

// ---------------------------------------


template <class Hamiltonian, class Params, class Callable>
void RegularLatticeLoop(RegistryDB& reg, const std::string outbasedir, const std::vector<Params>& hp, Callable filter)
{
	const auto name      = reg.Get<std::string>("mc", "General", "Hamiltonian" );
	      auto nreplicas = reg.Get<std::vector<int>>("mc", name, "rep" );
	const auto nL  	 = reg.Get<std::vector<int>>("mc", name, "L" );
	const auto dim 	 = reg.Get<int>("mc", name, "dim" );

	if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

	// lattice size loop
	for (std::size_t j=0; j<nL.size(); j++)
	{
		// prepare
		int L = nL[j];
		cout << endl << "L = " << L << endl << endl;

		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";

        	MARQOVConfig mp(outpath);
        	mp.setnsweeps(2);
		mp.setncluster(int(L/2));
		mp.setwarmupsteps(300);
		mp.setgameloopsteps(600);

		makeDir(mp.outpath);

		auto params = finalize_parameter_pair(mp, hp);
		auto rparams = replicator_pair(params, nreplicas[j]);

		// lattice
		RegularHypercubic latt(L, dim);

		// set up and execute
 		auto f = [&filter, &latt, &outbasedir, L](auto p){return filter(latt, p);}; //partially apply filter
 		Loop<Hamiltonian, RegularHypercubic>(rparams, f);
//        std::cout<<"End in RegularLatticeLoop "<<std::endl;
//        exit(0);
	}
}


// ---------------------------------------

std::string selectsim_startup(RegistryDB& registry)
{
	const auto ham        = registry.Get<std::string>("mc", "General", "Hamiltonian" );
	const auto dim 	  = registry.Get<int>("mc", ham, "dim" );
	const auto nreplicas  = registry.Get<std::vector<int>>("mc", ham, "rep" );
	const auto nreplicass = registry.Get<std::string>("mc", ham, "rep" );
	const auto nL  	  = registry.Get<std::vector<int>>("mc", ham, "L" );
	const auto nLs 	  = registry.Get<std::string>("mc", ham, "L" );

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
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);
 		RegularLatticeLoop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
	}
	else if (ham == "Heisenberg")
	{
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);
		RegularLatticeLoop<Heisenberg<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}
    else if (ham == "Phi4")
    {
		auto beta   = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto lambda = registry.Get<std::vector<double> >("mc", ham, "lambda");
		auto mass   = registry.Get<std::vector<double> >("mc", ham, "mass");

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
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto D    = registry.Get<std::vector<double> >("mc", ham, "D");
		auto parameters = cart_prod(beta, J, D);

		write_logfile(registry, beta);
 		RegularLatticeLoop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
    }
    else if (startswith(ham, "AshkinTeller"))
    {
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto K    = registry.Get<std::vector<double> >("mc", ham, "K");
		auto parameters = cart_prod(beta, J, K);

		write_logfile(registry, beta);
 		RegularLatticeLoop<AshkinTeller<int>>(registry, outbasedir, parameters, defaultfilter);
	}
	else if (ham == "XXZAntiferro")
	{
		auto beta     = registry.Get<std::vector<double>>("mc", ham, "beta");
		auto extfield = registry.Get<std::vector<double>>("mc", ham, "extfield");
		auto aniso    = registry.Get<std::vector<double>>("mc", ham, "aniso");
		auto parameters = cart_prod(beta, aniso, extfield);

		write_logfile(registry, beta);
		RegularLatticeLoop<XXZAntiferro<double, double> >(registry, outbasedir, parameters, defaultfilter);
	}
	else if (ham == "XXZAntiferroSingleAniso")
	{
		auto beta        = registry.Get<std::vector<double>>("mc", ham, "beta");
		auto extfield    = registry.Get<std::vector<double>>("mc", ham, "extfield");
		auto aniso       = registry.Get<std::vector<double>>("mc", ham, "aniso");
		auto singleaniso = registry.Get<std::vector<double>>("mc", ham, "singleaniso");
		auto parameters = cart_prod(beta, extfield, aniso, singleaniso);

		write_logfile(registry, extfield);
		RegularLatticeLoop<XXZAntiferroSingleAniso<double,double> >(registry, outbasedir, parameters, xxzfilter);
	}
	else if (startswith(ham, "Bimodal-Ising-EdwardsAnderson"))
	{
		const auto ham        = registry.Get<std::string>("mc", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc", ham, "L" );

		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);

	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
	        	MARQOVConfig mp(outpath);
	        	mp.setnsweeps(15);
			mp.setncluster(0);
			mp.setwarmupsteps(300);
			mp.setgameloopsteps(300);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

			// perform simulations
		 	Loop< EdwardsAndersonIsing<int>, RegularRandomBond<BimodalPDF>>(rparams, defaultfilter_triple);
		}
	}
	else if (startswith(ham, "Gaussian-Ising-EdwardsAnderson"))
	{
		const auto ham        = registry.Get<std::string>("mc", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc", ham, "L" );

		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);

	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
	        	MARQOVConfig mp(outpath);
	        	mp.setnsweeps(50);
			mp.setncluster(0);
			mp.setwarmupsteps(100);
			mp.setgameloopsteps(1000);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

			// perform simulations
		 	Loop< EdwardsAndersonIsing<int>, RegularRandomBond<GaussianPDF>>(rparams, defaultfilter_triple);
		}
	}
	/*
	else if (ham == "SSH")
	{

		auto beta   = registry.Get<std::vector<double> >("mc", ham, "betaMC");
		auto betaQM = registry.Get<double>("mc", ham, "betaQM");
		auto m      = registry.Get<std::vector<double> >("mc", ham, "m");
		auto k      = registry.Get<std::vector<double> >("mc", ham, "k");

		const auto name      = registry.Get<std::string>("mc", "General", "Hamiltonian" );
		      auto nreplicas = registry.Get<std::vector<int>>("mc", name, "rep" );
		const auto nL 	      = registry.Get<std::vector<int>>("mc", name, "L" );
		const auto nLtime  	 = registry.Get<std::vector<int>>("mc", name, "Ltime" );
		const auto dim 	 = registry.Get<int>("mc", name, "dim" );



		// set up replicas
		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			for (std::size_t jj=0; jj<nLtime.size(); jj++)
			{
				// prepare
				int L = nL[j];
				int Ltime  = nLtime[jj];
				cout << endl << "L_space = " << L << "\t" << "L_time = " << Ltime << endl << endl;
	
				std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
	
	     	   	MARQOVConfig mp(outpath);
	     	   	mp.setnsweeps(5);
				mp.setncluster(0);
				mp.setwarmupsteps(1000);
				mp.setgameloopsteps(5000);

				mp.outname = "Ltime"+std::to_string(Ltime);
	
				makeDir(mp.outpath);

	
				// compute delta tau
				std::vector<double> dtau = {betaQM/double(Ltime)};

				// set up parameters
				auto hp = cart_prod(beta, m, k, dtau);
				auto params = finalize_parameter_pair(mp, hp);
				auto rparams = replicator_pair(params, nreplicas[j]);
	
				// lattice
				SSHLattice latt(L, Ltime, dim);
	
				// set up and execute
	 			auto f = [&latt, &outbasedir, L](auto p){return sshfilter(latt, p);}; //partially apply filter
	 			Loop<SSH<double>, SSHLattice>(rparams, f);
			}
		}
	}
	*/
	else if (ham == "IsingCC")
	{
		const auto ham        = registry.Get<std::string>("mc", "General", "Hamiltonian" );
		const auto dim 	  = registry.Get<int>("mc", ham, "dim" );
		      auto nreplicas  = registry.Get<std::vector<int>>("mc", ham, "rep" );
		const auto nL  	  = registry.Get<std::vector<int>>("mc", ham, "L" );

		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }

		auto beta = registry.Get<std::vector<double> >("mc", "IsingCC", "beta");
		auto J    = registry.Get<std::vector<double> >("mc", "IsingCC", "J");

		auto hp = cart_prod(beta, J);
		write_logfile(registry, beta);

	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare output
			int L = nL[j];
			cout << endl << "L = " << L << endl << endl;
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
			makeDir(outpath);
	
			// Monte Carlo parameters
	        	MARQOVConfig mp(outpath);
	        	mp.setnsweeps(5);
			mp.setncluster(15);
			mp.setwarmupsteps(500);
			mp.setgameloopsteps(1500);

			// lattice parameters
			auto lp = std::make_tuple(L,dim);

			// form parameter triple and replicate
			auto params  = finalize_parameter_triple(lp, mp, hp);
			auto rparams = replicator(params, nreplicas[j]);

			// perform simulations
		 	Loop<Ising<int>, ConstantCoordinationLattice<Poissonian>>(rparams, defaultfilter_triple);
		}
	}
    else if (ham == "IrregularIsing1")
    {
		//
		// construct irregular lattice and pass it to MARQOV as a reference
		//

		const int L = 32;
		const int dim = 2;
		ConstantCoordinationLattice<Poissonian> ccl(L,dim);

		// prepare output
		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		makeDir(outpath);

		// Hamiltonian parameters
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		std::vector<double> myj = {-1.0};
		auto hp = cart_prod(beta, myj);

		// Monte Carlo parameters
		MARQOVConfig mp(outpath);
		mp.setrepid(1);
		mp.setnsweeps(L);
		mp.setncluster(10);

		auto params = finalize_parameter_pair(mp, hp);

		// partially apply filter
		auto f = [&ccl](auto p){return defaultfilter(ccl, p);};

		// perform simulations
		Loop<Ising<int>, ConstantCoordinationLattice<Poissonian>>(params, f);
	}
    else if (ham == "BlumeCapelBipartite")
    {
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto DA   = registry.Get<std::vector<double> >("mc", ham, "DA");
		auto DB   = registry.Get<std::vector<double> >("mc", ham, "DB");
		auto parameters = cart_prod(beta, J, DA, DB);

		const auto name      = registry.Get<std::string>("mc", "General", "Hamiltonian" );
		      auto nreplicas = registry.Get<std::vector<int>>("mc", name, "rep" );
		const auto nL 	      = registry.Get<std::vector<int>>("mc", name, "L" );
		const auto dim 	 = registry.Get<int>("mc", name, "dim" );

		write_logfile(registry, beta);

		// set up replicas
		if (nreplicas.size() == 1) { for (int i=0; i<nL.size()-1; i++) nreplicas.push_back(nreplicas[0]); }
	
		// lattice size loop
		for (std::size_t j=0; j<nL.size(); j++)
		{
			// prepare
			int L = nL[j];
	
			std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
	
	     	MARQOVConfig mp(outpath);
        		mp.setnsweeps(2);
			mp.setncluster(int(L/2));
			mp.setwarmupsteps(500);
			mp.setgameloopsteps(2500);

			makeDir(mp.outpath);

			// set up parameters
			auto params = finalize_parameter_pair(mp, parameters);
			auto rparams = replicator_pair(params, nreplicas[j]);

			// lattice
			SimpleBipartite latt(L, dim);



			// test area
//			auto terms = get_terms<SimpleBipartite>(latt, 0);
//			cout << "---> " << terms[0] << endl << endl;

//			RegularHypercubic latt2(L, dim);
//			auto terms2 = get_terms<RegularHypercubic>(latt2, 0);
//			cout << "---> " << terms2[0] << endl << endl;



	
			// set up and execute
	 		auto f = [&latt, &outbasedir, L](auto p){return defaultfilter(latt, p);}; //partially apply filter
	 		Loop<BlumeCapelBipartite<int>, SimpleBipartite>(rparams, f);
		}
	}
}







int main()
{

	// read config files
	RegistryDB registry("../src/config");

	// remove old output and prepare new one
	std::string outbasedir = registry.Get<std::string>("mc", "IO", "outdir" );
	std::string logbasedir = registry.Get<std::string>("mc", "IO", "logdir" );

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
