#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <tuple>
#include "rndwrapper.h"
#include "regular_lattice.h"
#include "vectorhelpers.h"
#include "helpers.h"
#include "cartprod.h"
#include "registry.h"


#include <iomanip>      // std::setprecision
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

#include "marqov.h"
#include "neighbourclass.h"
#include "svmath.h"

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"
#include "XXZAntiferroSingleAniso.h"
#include "AshkinTeller.h"

#include "systemtools.h"

template <class T1, class T2, class T3>
class Triple
{
public:
    Triple(T1 t1, T2 t2, T3 t3) : first(t1), second(t2), third(t3) {}
    T1 first;
    T2 second;
    T3 third;
};

template< class T1, class T2, class T3 >
constexpr auto make_triple( T1&& t, T2&& u, T3&& v) {return Triple<typename std::decay<T1>::type, typename std::decay<T2>::type, typename std::decay<T3>::type>(t,u,v);}

//two examples on how to extend the parsing capabilities of the registry
template <typename T>
struct GetTrait<std::vector<T> >//helper trait to break up a string at various predefined seperators
{
    static std::vector<T> Convert(std::string& arg)
    {
        const std::string delim("; ,");
        std::vector<T> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(GetTrait<T>::Convert(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) ));
        } while ((posold = arg.find_first_not_of(delim, pos) ) != std::string::npos);
        return retval;
    }
};

template <>
struct GetTrait<std::vector<std::string> >//helper trait to break up a string at various predefined seperators
{
    static std::vector<std::string> Convert(std::string arg)
    {
        const std::string delim("; ,");
        std::vector<std::string> retval;
        std::size_t pos = 0;
        std::size_t posold = 0;
        do
        {
            retval.push_back(arg.substr(posold, ( pos = arg.find_first_of(delim, posold) ) - posold) );
        } while ((posold = arg.find_first_not_of(delim, pos) ) != std::string::npos);
        return retval;
    }
};

void write_logfile(RegistryDB& reg, std::vector<double> loopvar)
{
	std::string logdir  = reg.Get<std::string>("mc", "IO", "logdir" );
	std::string logfile = reg.Get<std::string>("mc", "IO", "logfile" );
	std::ofstream os(logdir+"/"+logfile);
	os << std::setprecision(7);
	for (std::size_t i=0; i<loopvar.size(); i++) os << loopvar[i] << endl;
	os.close();
}

// //C++17 make_from_tuple from cppreference adapted for emplace.
// template <class Cont, class Latt, class Tuple, std::size_t... I>
// constexpr auto emplace_from_tuple_impl(Cont&& cont, Latt&& latt, MARQOV::MARQOVConfig&& mc, Tuple&& t, std::index_sequence<I...> )
// {
//   return cont.emplace_back(std::forward<Latt>(latt), std::forward<decltype(mc)>(mc), std::get<I>(std::forward<Tuple>(t))...) ;
// }
// 
// /** A function to construct an object in a container directly from a tuple
//  * @param cont the container where we append to.
//  * @param t the tuple containing the arguments.
//  */
// template <class Cont, class T, class Tuple, class Latt>
// constexpr auto emplace_from_tuple(Cont&& cont, Latt&& latt, T&& mc, Tuple&& t )
// {
//     return emplace_from_tuple_impl(cont, std::forward<Latt>(latt), std::forward<T>(mc), std::forward<Tuple>(t),
//         std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
// }

//c++17 make_from_tuple from cppreference adapted for emplace
template <class Cont, class Tuple1, class Tuple2, std::size_t... I>
constexpr auto emplace_from_tuple_impl(Cont&& cont, Tuple1&& t1, MARQOV::MARQOVConfig&& mc, Tuple2&& t2, std::index_sequence<I...> )
{
  return cont.emplace_back(std::forward<Tuple1>(t1), std::forward<MARQOV::MARQOVConfig>(mc), std::get<I>(std::forward<Tuple2>(t2))...) ;
}
 
template <class Cont, class Tuple1, class Tuple2>
constexpr auto emplace_from_tuple(Cont&& cont, Tuple1&& t1, MARQOV::MARQOVConfig&& mc, Tuple2&& t2 )
{
    return emplace_from_tuple_impl(cont, std::forward<Tuple1>(t1), std::forward<MARQOV::MARQOVConfig>(mc), std::forward<Tuple2>(t2),
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
        emplace_from_tuple(	sims, 
	   					std::forward<decltype(std::get<0>(t1))>(std::get<0>(t1)), 
						std::forward<MARQOV::MARQOVConfig>(std::get<1>(t1)), std::get<2>(t1));
    }
}

// ---------------------------------------


using namespace MARQOV;


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
auto createsims(std::vector<Triple<LArgs, MARQOVConfig, HArgs> >& params, Callable c)
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
void loop(MARQOVConfig& mc, const std::vector<Parameters>& hamparams, Callable filter)
{
	// number of EMCS during relaxation and measurement
	mc.setwarmupsteps(1500);
	mc.setgameloopsteps(5000);

	std::vector<std::pair<MARQOVConfig, Parameters> > params;
	for(std::size_t i = 0; i < hamparams.size(); ++i)
	{
	    auto mc2(mc);
	    mc2.setid(i);
	    params.push_back(make_pair(mc2, hamparams[i]));
	}

	auto sims = createsims<Hamiltonian, Lattice>(params, filter);

	// perform simulation
	#pragma omp parallel for
	for(std::size_t i = 0; i < sims.size(); ++i)
	{
		auto& marqov = sims[i];
		marqov.init();
//		marqov.init_hot();
//		marqov.init_cold_Ising_like(); cout << "AT init is on!!!" << endl;
		marqov.wrmploop();
		marqov.gameloop();
	}
}


// ---------------------------------------


template <class Hamiltonian, class Params, class Callable>
void RegularLatticeloop(RegistryDB& reg, const std::string outbasedir, const std::vector<Params>& parameters, Callable filter)
{
	const auto dim 		= reg.Get<int>("mc", "General", "dim" );
	const auto nreplicas 	= reg.Get<int>("mc", "General", "nreplicas" );
	const auto nL  		= reg.Get<std::vector<int>>("mc", "General", "nL" );
	const std::string name 	= reg.Get<std::string>("mc", "General", "Hamiltonian" );
	const std::string nLs 	= reg.Get<std::string>("mc", "General", "nL" );

	cout << endl;
	cout << "Hamiltonian: \t" << name << endl;
	cout << "Dimension: \t" << dim << endl;
	cout << "Lattice sizes:\t" << nLs << endl;
	cout << "Replicas:\t" << nreplicas << endl;


	// lattice size loop
	for (std::size_t j=0; j<nL.size(); j++)
	{
		// prepare
		int L = nL[j];
		cout << endl << "L = " << L << endl << endl;

		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";

        	MARQOVConfig mc(outpath);
        	mc.setnsweeps(5);
		mc.setncluster(15);

		makeDir(mc.outpath);

		// lattice
		RegularLattice latt(L, dim);

		// set up and execute
 		auto f = [&filter, &latt, &outbasedir, L](auto p){return filter(latt, p);}; //partially apply filter
 		loop<Hamiltonian, RegularLattice>(mc, parameters, f);
	}
}


// ---------------------------------------


void selectsim(RegistryDB& registry, std::string outbasedir, std::string logbasedir)
{
	// extract Hamiltonian type and number of replicas
	const std::string ham = registry.Get<std::string>("mc", "General", "Hamiltonian" );
	const int   nreplicas = registry.Get<int>("mc", "General", "nreplicas" );

	// by construction temperatures have to go first in cart_prod, but sometimes
	// we want it sorted by "id", therefore the two are swapped afterwards
	std::vector<double> id(nreplicas);
	for (int i=0; i<nreplicas; ++i) id[i] = i;

	std::vector<MARQOVConfig> mcs(nreplicas, MARQOVConfig(outbasedir));
	for (int i = 0; i < nreplicas; ++i) mcs[i].setid(i);



	// filter to determine output file path and name
	// the filter _must_ set p.first.outname!
	auto defaultfilter = [](auto& latt, auto p)
	{
		auto& mp = p.first;		// Monte Carlo params
		auto& hp = p.second;	// Hamiltonian params
		
		std::string str_id    = std::to_string(mp.id);
		std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
		mp.outname = str_beta;
		return std::tuple_cat(std::forward_as_tuple(latt), p);
	};



	if (ham == "Ising")
	{
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto parameters = cart_prod(beta, J);
		//FIXME: is this really intended to have nreplicas of each (beta, J) -> yes ;)

		write_logfile(registry, beta);
 		RegularLatticeloop<Ising<int>>(registry, outbasedir, parameters, defaultfilter);
	}
    else if (ham == "Heisenberg")
    {
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);
		RegularLatticeloop<Heisenberg<double, double> >(registry, outbasedir, parameters, defaultfilter);
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
		RegularLatticeloop<Phi4<double, double> >(registry, outbasedir, parameters, defaultfilter);
    }
    else if (ham == "BlumeCapel")
    {
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto D    = registry.Get<std::vector<double> >("mc", ham, "D");
		auto parameters = cart_prod(beta, J, D);

		write_logfile(registry, beta);
 		RegularLatticeloop<BlumeCapel<int>>(registry, outbasedir, parameters, defaultfilter);
    }
    else if (ham == "AshkinTeller")
    {
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto K    = registry.Get<std::vector<double> >("mc", ham, "K");
		auto parameters = cart_prod(beta, J, K);

		write_logfile(registry, beta);
 		RegularLatticeloop<AshkinTeller<int>>(registry, outbasedir, parameters, defaultfilter);
    }
    else if(ham == "XXZAntiferro")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj, myj, myj);
        RegularLatticeloop<XXZAntiferro<double, double> >(registry, outbasedir, parameters, defaultfilter);
    }
    else if(ham == "XXZAntiferroSingleAniso")
    {
		auto beta        = registry.Get<std::vector<double>>("mc", ham, "beta");
		auto extfield    = registry.Get<std::vector<double>>("mc", ham, "extfield");
		auto aniso       = registry.Get<std::vector<double>>("mc", ham, "aniso");
		auto singleaniso = registry.Get<std::vector<double>>("mc", ham, "singleaniso");

		auto parameters = cart_prod(id, beta, extfield, aniso, singleaniso);

		for (auto& param_tuple : parameters) // swap "id" and "temperature"
					std::swap(std::get<0>(param_tuple), std::get<1>(param_tuple));

		write_logfile(registry, extfield);


		// write a filter to determine output file path and name
		auto xxzfilter = [](RegularLattice& latt, auto p)
		{	
			auto& mp = p.first;		// Monte Carlo params
			auto& hp = p.second;	// Hamiltonian params
		
			std::string str_id    = std::to_string(int(std::get<1>(hp)));
			std::string str_beta  = "beta"+std::to_string(std::get<0>(hp));
			std::string str_extf  = "extf"+std::to_string(std::get<2>(hp));
			
			mp.outname = str_beta+"_"+str_extf+"_"+str_id;

			return std::tuple_cat(std::forward_as_tuple(latt), p);
		};

		RegularLatticeloop<XXZAntiferroSingleAniso<double,double> >(registry, outbasedir, parameters, xxzfilter);
	}
    else if(ham == "IrregularIsing")
    {
		std::vector<std::vector<int>> dummy;
		Neighbours<int32_t> nbrs(dummy);

		auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
		std::vector<double> myj = {1.0};
		auto hamparams = cart_prod(betas, myj);

		// extract lattice size and prepare directories
		int L = nbrs.size();
		cout << endl << "L = " << L << endl << endl;

		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		MARQOVConfig mc(outpath);
		makeDir(mc.outpath);

		mc.setid(1);
		mc.setnsweeps(2*L);
		mc.setncluster(1);

		// lattice
		auto f = [&defaultfilter, &nbrs](auto p){return defaultfilter(nbrs, p);}; //partially apply filter
		loop<Ising<int>, Neighbours<int32_t> >(mc, hamparams, f);
	}
	else if(ham == "IrregularIsing2")
	{
		auto beta = registry.Get<std::vector<double> >("mc", ham, "beta");
		auto J    = registry.Get<std::vector<double> >("mc", ham, "J");
		auto parameters = cart_prod(beta, J);

		write_logfile(registry, beta);

		auto otherfilter = [](auto p)
		{
			// write a filter to determine output file path and name
            	auto& lp = p.first;
            	auto& mp = p.second;
            	auto& hp = p.third;

 			auto str_id    = std::to_string(mp.id);
			auto str_beta  = "beta"+std::to_string(std::get<0>(hp));
			auto str_L     = std::to_string(std::get<0>(lp));

			mp.outname = str_beta+"_"+str_id;

			return p;
		};


		const int L   = 42;
		const int dim = 2;
		
		std::string outpath = outbasedir+"/"+std::to_string(L)+"/";
		MARQOVConfig mc(outpath);
		makeDir(mc.outpath);

		auto t = make_triple(std::make_tuple(L,dim), mc, parameters[0]);
		std::vector<decltype(t)> p = {t};

		auto sims = createsims<Ising<int>, RegularLattice >(p, otherfilter);

		// perform simulation
		#pragma omp parallel for
		for(std::size_t i = 0; i < sims.size(); ++i)
		{
			auto& marqov = sims[i];
	
			marqov.init();
//			marqov.init_hot();
			marqov.wrmploop();
			marqov.gameloop();
		}

//        auto t = make_pair(std::make_tuple(dummy), std::tuple_cat(std::make_tuple(outdir), parameters[0]));
//        createsims<Ising<int>, Neighbours<int32_t> >(p, otherfilter);
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
