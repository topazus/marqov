#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <cstdlib>
#include <fstream>
#include "rndwrapper.h"
#include "regular_lattice.h"
#include "vectorhelpers.h"
#include "registry.h"

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

		// temperature loop
		for (int i=0; i<nbeta; i++)
		{
			double currentbeta = betastart + i*betastep; 
			cout << "beta = " << currentbeta << endl;

			// set up lattice
			RegularLattice lattice(L, dim);
		
			// set up outfile
			std::string outfile = outdir+"/"+std::to_string(L)+"/"+std::to_string(i)+".h5";

			// set up model
            {Marqov<RegularLattice, Ising<int> > marqov(lattice, outfile, currentbeta, -1.0);}
			{Marqov<RegularLattice, BlumeCapel<int> > marqov(lattice, outfile, currentbeta, currentbeta);}
			{Marqov<RegularLattice, Heisenberg<double,double> > marqov(lattice, outfile, currentbeta, -1.0);}
			{Marqov<RegularLattice, Phi4<double,double> > marqov(lattice, outfile, currentbeta, currentbeta, -4.5);}
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
