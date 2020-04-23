#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "rndwrapper.h"
#include "regular_lattice.h"
#include "vectorhelpers.h"
#include "helpers.h"
#include "cartprod.h"
#include "registry.h"


//helper, delete later
#include <typeinfo>
#include <cxxabi.h>

#include <iomanip>      // std::setprecision
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

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


// ---------------------------------------

const int myid = 0; // remove once a parallelization is available

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"
#include "XXZAntiferroSingleAniso.h"
#include "AshkinTeller.h"

#include "marqov.h"
#include "disordered_lattice.h"

int main()
{
	// read config files
	RegistryDB registry("../src/config");


	// remove old output and prepare new one
	std::string outdir = registry.Get<std::string>("mc", "IO", "outdir" );
	std::string logdir = registry.Get<std::string>("mc", "IO", "logdir" );

	std::string command;
	command = "rm -r " + outdir;
	system(command.c_str());
	command = "rm -r " + logdir;
	system(command.c_str());

	makeDir(outdir);
	makeDir(logdir);


	// ------------------ live view -----------------------
	/*

		const int L = 30;
		RegularLattice lattice(L, 3);
		std::string outfile = outdir+"temp.h5";
		Marqov<RegularLattice, XXZAntiferro<double,double> > marqov(lattice, outfile, 1/0.36);
		marqov.init_hot();
		const int ncluster = 0;
		const int nsweeps  = L/2; 
		marqov.wrmploop(50, ncluster, nsweeps);
		marqov.gameloop_liveview();

	*/
	// ----------------------------------------------------



	// --------- unpack configuration file ---------

	const auto dim       = registry.Get<int>("mc", "General", "dim" );
	const auto nL        = registry.Get<std::vector<int> >("mc", "General", "nL" );
	const int  nreplicas = registry.Get<int>("mc", "General", "nreplicas" );

	auto loopstyle = registry.Get<std::string>("mc", "General", "loopstyle" );

	double lvstart = registry.Get<double>("mc", "General", "lvstart" );
	double lvfinal = registry.Get<double>("mc", "General", "lvfinal" );
	int    lvsteps = registry.Get<int>("mc", "General", "lvsteps" );


	// create range for loop variable
	std::vector<double> loopvar = create_range(lvstart, lvfinal, lvsteps, "linear", 0);

	// write values in external fields in logfile
	std::string logfile = registry.Get<std::string>("mc", "IO", "logfile" );
	std::ofstream os(logdir+"/"+logfile);
		os << std::setprecision(7);
		for (int i=0; i<loopvar.size(); i++) os << loopvar[i] << endl;
	os.close();



	// ---------------- main loop -----------------

	cout << endl << "The dimension is " << dim << endl;

	// lattice size loop
	for (int j=0; j<nL.size(); j++)
	{
		// extract lattice size and prepare directories
		const int L = nL[j];
		cout << endl << "L = " << L << endl << endl;
		makeDir(outdir+"/"+std::to_string(L));

		for (int k=0; k<lvsteps; k++)
		{

//			const double beta = 1/1.5;
//			const double p = loopvar[k];
			
			const double beta = loopvar[k];

        
			#pragma omp parallel for
			for(std::size_t i = 0; i < nreplicas; ++i)
			{
				std::string str_var = "var"+std::to_string(beta);
				std::string str_id   = std::to_string(i);
				std::string outname   = str_var+"_"+str_id+".h5";
				std::string outsubdir = outdir+"/"+std::to_string(L)+"/";


				// lattice and model
//				RegularRandomBond<double> latt(dim, L, p);
//				Marqov<RegularRandomBond<double>, Ising<int>> sim(latt, outsubdir+outname, beta);

				RegularLattice latt(L, dim);
				Marqov<RegularLattice, AshkinTeller<int>> sim(latt, outsubdir+outname, beta);


				// number of cluster updates and metropolis sweeps
				const int ncluster = 1;
				const int nsweeps  = 10;
				
				// number of EMCS during relaxation and measurement
				const int nrlx = 500;
				const int nmsr = 1500;  
				
				
				// perform simulation
//				sim.init_hot();
				sim.init_cold_Ising_like();
				sim.wrmploop(nrlx, ncluster, nsweeps, i);
				sim.gameloop(nmsr, ncluster, nsweeps, i);
			}
		}
	}
}
