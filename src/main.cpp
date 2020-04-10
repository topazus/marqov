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
#include "cartprod.h"
#include "registry.h"

//helper, delete later
#include <typeinfo>
#include <cxxabi.h>

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
    {
        auto temp = c(p);
        emplace_from_tuple(sims, temp);
    }
}

// ---------------------------------------

const int myid = 0; // remove once a parallelization is available

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"

#include "marqov.h"

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


	// extract Monte Carlo parameters
	auto  nL = registry.Get<std::vector<int> >("mc", "General", "nL" );
	int    nbeta     = registry.Get<int>("mc", "General", "nbeta" );
	double betastart = registry.Get<double>("mc", "General", "betastart" );
	double betaend   = registry.Get<double>("mc", "General", "betaend" );
	double betastep = (betaend - betastart) / double(nbeta);


	// write temperatures in logfile
	std::string logfile = registry.Get<std::string>("mc", "IO", "logfile" );
	std::ofstream os(logdir+"/"+logfile);
	for (int i=0; i<nbeta; i++) os << betastart + i*betastep << endl;
	os.close();



	// lattice dimension
	const int dim = 3;

	// lattice size loop
	for (int j=0; j<nL.size(); j++)
	{
		// extract lattice size and prepare directories
		const int L = nL[j];
		cout << endl << "L = " << L << endl << endl;
		cout << outdir+std::to_string(L) << endl;
		makeDir(outdir+std::to_string(L));
        
        nbeta=7;
        betastep = 0.1;
        // let's create some parameter vectors
        std::vector<double> anisos = {0.8, 1.0, 1.2};
        std::vector<double> anisos2 = {10,20,30,40,50};
        std::vector<double> anisos3 = {100,200,300,400,500};
        std::vector<double> betas(nbeta);
        for (int i = 0; i < nbeta; ++i)
            betas[i] = betastart + i*betastep;
        std::cout<<betas.size()<<" "<<betas[0]<<" "<<betas.back()<<std::endl;
        auto parameters = cart_prod(betas, anisos);
        std::vector<Marqov<RegularLattice, XXZAntiferro<double,double> >> sims;
        sims.reserve(parameters.size());//MARQOV has issues with copying
        RegularLattice latt(L, dim);
        fillsims(parameters, sims, 
                 [&latt, &outdir, L]( decltype(parameters[0]) p) {
                     std::string outfile = outdir+std::to_string(L)+"/beta"+std::to_string(std::get<0>(p))
                     +"aniso"+std::to_string(std::get<1>(p))+".h5";
                     std::cout<<outfile<<std::endl;
                     return std::tuple_cat(std::make_tuple(latt, outfile), p);
                });
        std::cout<< sims.size()<<std::endl;
        exit(-1);
		// temperature loop
		for (int i=0; i<nbeta; i++)
		{
			double currentbeta = betastart + i*betastep; 
			cout << "beta = " << currentbeta << endl;

			// set up lattice
			RegularLattice lattice(L, dim);
		
			// set up outfile
			std::string outfile = outdir+std::to_string(L)+"/"+std::to_string(i)+".h5";

			// set up model
//             {Marqov<RegularLattice, Ising<int> > marqov(lattice, outfile, currentbeta, -1.0);}
// 			{Marqov<RegularLattice, BlumeCapel<int> > marqov(lattice, outfile, currentbeta, currentbeta);}
// 			{Marqov<RegularLattice, Heisenberg<double,double> > marqov(lattice, outfile, currentbeta, -1.0);}
// 			{Marqov<RegularLattice, Phi4<double,double> > marqov(lattice, outfile, currentbeta, currentbeta, -4.5);}
			Marqov<RegularLattice, XXZAntiferro<double,double> > marqov(lattice, outfile, currentbeta, 1.0);


			// number of cluster updates and metropolis sweeps
			const int ncluster = 0;
			const int nsweeps  = L; 


			// number of EMCS during relaxation and measurement
			const int nrlx = 2500;
			const int nmsr = 5000;  


			// perform simulation
			marqov.init_hot();
			marqov.wrmploop(nrlx, ncluster, nsweeps);
			marqov.gameloop(nmsr, ncluster, nsweeps);

		}
	}
}
