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

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"
#include "XXZAntiferroSingleAniso.h"
#include "AshkinTeller.h"

#include "marqov.h"
#include "neighbourclass.h"


//helper, delete later
#include <typeinfo>
#include <cxxabi.h>


#include "systemtools.h"

//two examples on how to extend the parsing capapbilities of the registry
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

class Grid 
{
    std::vector<int> getnbr(int i);
};

// ------- elementary state vector calculus

template <class StateVector>
StateVector operator + (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] += rhs[i];
    return res;
}

template <class StateVector>
StateVector operator - (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] -= rhs[i];
    return res;
}

inline double mult(const double& a, const double& b)
{
    return a*b;
}

inline double mult(const double& a, const int& b)
{
    return a*double(b);
}

template <class VecType, class StateVector>
inline StateVector mult(const VecType& a, const StateVector& b)
{
    StateVector retval(b);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    retval[i] *= a[i];
    return retval;
}


inline double dot(const double& a, const double& b)
{
    return a*b;
}

template<class VecType>
inline typename VecType::value_type dot(const VecType& a, const VecType& b)
{
    typedef typename VecType::value_type FPType;
    return std::inner_product(begin(a), end(a), begin(b), 0.0);
}


template <class StateVector>
inline void reflect(StateVector& vec, const StateVector mirror)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	const double dotp = dot(vec,mirror);

	for (int i=0; i<SymD; i++) vec[i] -= 2*dotp*mirror[i];
}	

template <class Container>
inline void normalize(Container& a)
{
	typename Container::value_type tmp_abs=std::sqrt(dot(a, a));

	for (int i = 0; i < a.size(); ++i) a[i] /= tmp_abs;
}



template <class StateVector>
inline void coutsv(StateVector& vec)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	for (int i=0; i<SymD; i++) cout << vec[i] << "\t";
	cout << endl;
}

//c++17 make_from_tuple from cppreference adapted for emplace
template <class Cont, class Tuple, std::size_t... I>
constexpr auto emplace_from_tuple_impl(Cont&& cont, Tuple&& t, std::index_sequence<I...> )
{
  return cont.emplace_back(std::get<I>(std::forward<Tuple>(t))...) ;
}
 
template <class Cont, class Tuple>
constexpr auto emplace_from_tuple(Cont&& cont,  Tuple&& t )
{
    return emplace_from_tuple_impl(cont, std::forward<Tuple>(t),
        std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}

/** ToDo: introduce a filter..., e.g. a lambda*/
template <class Args, class T, class Grid>
void fillsims(const std::vector<Args>& args, std::vector<T>& sims, Grid& grid, std::string& fn)
{
    for(auto p : args)
    {
        auto args = std::tuple_cat(std::make_tuple(grid, fn), p);
        emplace_from_tuple(sims, args);
    }
}

template <class Args, class T, class Callable>
void fillsims(const std::vector<Args>& args, std::vector<T>& sims, Callable c)
{
    for(auto p : args)
        emplace_from_tuple(sims, c(p));
}

using namespace MARQOV;

template <class Hamiltonian, class Lattice, class Parameters, class Callable>
void loop(const std::vector<Parameters>& params, Callable filter, int nsweeps)
{
		// model
		std::vector<Marqov<Lattice, Hamiltonian> > sims;

		// simulation vector
		sims.reserve(params.size());//MARQOV has issues with copying -> requires reserve in vector

		fillsims(	params, sims, filter);

//init
        #pragma omp parallel for
		for(std::size_t i = 0; i < sims.size(); ++i) //for OMP
		{
			auto myid = i;
			auto& marqov = sims[i];

			// number of cluster updates and metropolis sweeps
			const int ncluster = 0;
			
			// number of EMCS during relaxation and measurement
			const int nrlx = 5000;
			// perform simulation
			marqov.init_hot();
			marqov.wrmploop(nrlx, ncluster, nsweeps, myid);
		}

std::cout<<"warmup done"<<std::endl;
		// ------------- execute -------------

		#pragma omp parallel for
		for(std::size_t i = 0; i < sims.size(); ++i) //for OMP
		{
			auto myid = i;
			auto& marqov = sims[i];

			// number of cluster updates and metropolis sweeps
			const int ncluster = 0;
			
			// number of EMCS during relaxation and measurement
			const int nmsr = 15000;  
			// perform simulation
			marqov.gameloop(nmsr, ncluster, nsweeps, myid);
		}
}

// ---------------------------------------

const int myid = 0; // remove once a parallelization is available


template <class Hamiltonian, class Params, class Callable>
void RegularLatticeloop(RegistryDB& reg, const std::string outdir, std::string logdir, const std::vector<Params>& parameters, Callable filter)
{
    const auto dim       = reg.Get<int>("mc", "General", "dim" );
	const auto nL        = reg.Get<std::vector<int> >("mc", "General", "nL" );

	// ---------------- main loop -----------------

	cout << endl << "The dimension is " << dim << endl;

	// lattice size loop
	for (int j=0; j<nL.size(); j++)
	{
		// extract lattice size and prepare directories
		int L = nL[j];
		cout << endl << "L = " << L << endl << endl;
		makeDir(outdir+"/"+std::to_string(L));
		// lattice
		RegularLattice latt(L, dim);
        auto f = [&filter, &latt, &outdir, L](auto p){return filter(latt, outdir, L, p);};//partially apply filter
        loop<Hamiltonian, RegularLattice>(parameters, f, 2*L); 
	}
}

void selectsim(RegistryDB& registry, std::string outdir, std::string logdir)
{
    const std::string ham = registry.Get<std::string>("mc", "General", "Hamiltonian" );
    const int  nreplicas = registry.Get<int>("mc", "General", "nreplicas" );
    auto defaultfilter = [](auto& latt, std::string outdir, int L, auto p)
    {
        // write a filter to determine output file path and name
        std::string str_beta = "beta"+std::to_string(std::get<0>(p));
        std::string outname   = str_beta+".h5";
        std::string outsubdir = outdir+"/"+std::to_string(L)+"/";
        return std::tuple_cat(std::forward_as_tuple(latt), std::make_tuple(outsubdir+outname), p);
    };
    if (ham == "Ising")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas); //, myj);
        RegularLatticeloop<Ising<int> >(registry, outdir, logdir, parameters, defaultfilter);
    }
    /*
    else if (ham == "Heisenberg")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj);
        RegularLatticeloop<Heisenberg<double, double> >(registry, outdir, logdir, parameters, defaultfilter);
    }
    else if (ham == "Phi4")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj, myj);
        RegularLatticeloop<Phi4<double, double> >(registry, outdir, logdir, parameters, defaultfilter);
    }
    else if (ham == "BlumeCapel")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj);
        RegularLatticeloop<BlumeCapel<int> >(registry, outdir, logdir, parameters, defaultfilter);
    }
    else if(ham == "XXZAntiferro")
    {
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj, myj,myj);
        RegularLatticeloop<XXZAntiferro<double, double> >(registry, outdir, logdir, parameters, defaultfilter);
    }
    else if(ham == "XXZAntiferroSingleAniso")
    {
        	// --------- unpack configuration file ---------

	auto loopstyle = registry.Get<std::string>("mc", "General", "loopstyle" );

	double lvstart = registry.Get<double>("mc", "General", "lvstart" );
	double lvfinal = registry.Get<double>("mc", "General", "lvfinal" );
	int    lvsteps = registry.Get<int>("mc", "General", "lvsteps" );

    auto parnames = registry.GetBlock("mc", ham).GetKeys();

	std::vector<std::vector<double>> par;
	for (int i=0; i<parnames.size(); i++)
	{
		auto parname = parnames[i];
		if (parname != lvname)
		{
			std::vector<double> parval = {0};
			parval[0] = registry.Get<double>("mc", ham, parnames[i]);
			par.push_back(parval);
		}
	}

	// create range for loop variable
	std::vector<double> loopvar = create_range(lvstart, lvfinal, lvsteps);
        // ------------ create parameter vector ---------------

		// todo: improve this section
		// by construction temperatures have to go first, but here 
		// we want it sorted by "id", therefore the two are swapped afterwards

		// replicas index (used as a fake Hamiltonian parameter)
		std::vector<double> id(nreplicas);
		for (int i=0; i<nreplicas; ++i) id[i] = i;

		// beta is loopvar
		auto parameters = cart_prod(id, loopvar, par[0], par[1], par[2]);

		// beta is not loopvar
//		auto beta = registry.Get<double>("mc", "General", "beta");
//		auto parameters = cart_prod(id, beta, loopvar, par[0], par[1]);
	

		// swap "id" and "temperature"
		for (auto& param_tuple : parameters) 
			std::swap(std::get<0>(param_tuple), std::get<1>(param_tuple));


	// write values in external fields in logfile
	std::string logfile = registry.Get<std::string>("mc", "IO", "logfile" );
	std::ofstream os(logdir+"/"+logfile);
		os << std::setprecision(7);
		for (int i=0; i<loopvar.size(); i++) os << loopvar[i] << endl;
	os.close();
                auto xxzfilter = [](RegularLattice& latt, std::string outdir, int L, decltype(parameters[0]) p)
				{
					// write a filter to determine output file path and name
					std::string str_beta = "beta"+std::to_string(std::get<0>(p));
					std::string str_extf = "extf"+std::to_string(std::get<2>(p));
					std::string str_id   = std::to_string(int(std::get<1>(p)));

					std::string outname   = str_beta+"_"+str_extf+"_"+str_id+".h5";
					std::string outsubdir = outdir+"/"+std::to_string(L)+"/";

					return std::tuple_cat(std::forward_as_tuple(latt), std::make_tuple(outsubdir+outname), p);
				};
            RegularLatticeloop<XXZAntiferroSingleAniso<double,double> >(registry, outdir, logdir, parameters, xxzfilter);
    }
    else if(ham == "IregularIsing")
    {
        std::vector<std::vector<int> > dummy;
        Neighbours<int32_t> nbrs(dummy);
        auto betas = registry.Get<std::vector<double> >("mc", ham, "betas");
        std::vector<double> myj = {1.0};
        auto parameters = cart_prod(betas, myj);
        		// extract lattice size and prepare directories
		int L = nbrs.size();
		cout << endl << "L = " << L << endl << endl;
		makeDir(outdir+"/"+std::to_string(L));
		// lattice
        auto f = [&defaultfilter, &nbrs, &outdir, L](auto p){return defaultfilter(nbrs, outdir, L, p);};//partially apply filter
        loop<Ising<int>, Neighbours<int32_t> >(parameters, f, 2*L);
    }
*/
}

int main()
{
	// read config files
	RegistryDB registry("../src/config");

	// remove old output and prepare new one
	std::string outdir = registry.Get<std::string>("mc", "IO", "outdir" );
	std::string logdir = registry.Get<std::string>("mc", "IO", "logdir" );

    //FIXME: NEVER DELETE USER DATA
	std::string command;
	command = "rm -r " + outdir;
	system(command.c_str());
	command = "rm -r " + logdir;
	system(command.c_str());

	makeDir(outdir);
	makeDir(logdir);
	
selectsim(registry, outdir, logdir);
}
