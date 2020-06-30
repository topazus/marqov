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
#include "geom/regular_lattice.h"
#include "geom/grid.h"
#include "geom/neighbourclass.h"
#include "geom/io.h"
#include "vectorhelpers.h"
#include "helpers.h"
#include "cartprod.h"
#include "registry.h"
#include "systemtools.h"
#include "replicate.h"
#include "marqov.h"
#include "svmath.h"

// Hamiltonians
#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"
#include "XXZAntiferroSingleAniso.h"
#include "AshkinTeller.h"

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


template<class ... Ts>
struct sims_helper {};

template <class H,  class L, class HArgstuple, size_t... S>
struct sims_helper<H, L, HArgstuple, std::index_sequence<S...> >
{
	typedef decltype(makeMarqov<H>(std::declval<L>(),
	                       		 std::declval<MARQOVConfig>(),
	                       		 std::declval<typename std::tuple_element<S, HArgstuple>::type>()...
							 )) MarqovType;
};

template <class H, class L, class HArgs, class LArgs>
struct sims_helper<H, L, Triple<LArgs, MARQOVConfig, HArgs> >
{
    typedef decltype(makeMarqov<H,L>(std::declval<MARQOVConfig>(),  
    							  std::declval<std::pair<LArgs, HArgs>& >()
							  )) MarqovType;
};


// ---------------------------------------


/** The case where Marqov allocates a lattice
 */
template <class H, class L, class LArgs, class HArgs, class Callable>
auto createsims(const std::vector<Triple<LArgs, MARQOVConfig, HArgs> >& params, Callable c)
{
    typedef typename sims_helper<H, L, Triple<LArgs, MARQOVConfig, HArgs> >::MarqovType MarqovType;  
    //create simulations
    std::vector<MarqovType> sims;
    sims.reserve(params.size());
    fillsims(params, sims, c);
    return sims;
}

/** The old case where the lattice is a reference passed in through the filter....
 * @param params parameters
 * @param c A filter
 */
template <class H, class L, class Args, class Callable>
auto createsims(const std::vector<std::pair<MARQOVConfig, Args>>& params, Callable c)
{
    std::size_t constexpr tsize = std::tuple_size<typename std::remove_reference<Args>::type>::value;
    typedef typename sims_helper<H, L, Args, std::make_index_sequence<tsize> >::MarqovType MarqovType;

    //create simulations
    std::vector<MarqovType> sims;
    sims.reserve(params.size());
    fillsims(params, sims, c);
    return sims;
}


// ---------------------------------------


template <class Hamiltonian, class Lattice, class Parameters, class Callable>
void Loop(const std::vector<Parameters>& params, Callable filter)
{
	auto sims = createsims<Hamiltonian, Lattice>(params, filter);

	// perform simulation
	#pragma omp parallel for
	for(std::size_t i = 0; i < sims.size(); ++i)
	{
		auto& marqov = sims[i];
		marqov.init();
//		marqov.debugloop(100,0,1);
		marqov.wrmploop();
		marqov.gameloop();
	}
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
        	mp.setnsweeps(5);
		mp.setncluster(15);
		mp.setwarmupsteps(500);
		mp.setgameloopsteps(1500);

		makeDir(mp.outpath);

		auto params = finalize_parameter_pair(mp, hp);
		auto rparams = replicator_pair(params, nreplicas[j]);

		// lattice
		RegularHypercubic latt(L, dim);

		// set up and execute
 		auto f = [&filter, &latt, &outbasedir, L](auto p){return filter(latt, p);}; //partially apply filter
 		Loop<Hamiltonian, RegularHypercubic>(rparams, f);
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


	// -------------------- filters --------------------

	// filter to determine output file path and name
	// the filter _must_ set the outname

	auto defaultfilter = [](auto& latt, auto p)
	{
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
		
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		mp.outname = str_beta+"_"+str_repid;
		return std::tuple_cat(std::forward_as_tuple(latt), p);
	};

	auto defaultfilter_triple = [](auto p)
	{
       	auto& lp = p.first;
       	auto& mp = p.second;
       	auto& hp = p.third;

 		auto str_repid = std::to_string(mp.repid);
		auto str_beta  = "beta"+std::to_string(std::get<0>(hp));
		auto str_L     = std::to_string(std::get<0>(lp));

		mp.outname = str_beta+"_"+str_repid;

		return p;
	};

	auto xxzfilter = [](auto& latt, auto p)
	{	
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
	
		std::string str_repid = std::to_string(mp.repid);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		std::string str_extf  = "extf"+std::to_string(std::get<1>(hp));
		mp.outname = str_beta+"_"+str_extf+"_"+str_repid;

		return std::tuple_cat(std::forward_as_tuple(latt), p);
	};





	// ----------------- select simulation ------------------

	if (ham == "Ising")
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
	/*
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
		auto f = [&defaultfilter, &ccl](auto p){return defaultfilter(ccl, p);};

		// perform simulations
		Loop<Ising<int>, ConstantCoordinationLattice<Poissonian>>(params, f);
	}
	*/
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
