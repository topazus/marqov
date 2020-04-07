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


class Grid 
{
    std::vector<int> getnbr(int i);
};



template <class StateVector>
class Interaction
{
	public:
	    double J;
	    virtual StateVector operator() (StateVector& phi_i) = 0;
};


template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
	    CouplingType h;
	    virtual CouplingType operator() (StateVector& phi) = 0;
	private:
};



template <class StateSpace, class StateVector>
class MultiSite
{
	public: 
	    const double k;
	    virtual double operator() (StateVector& sv, int svpos, StateSpace s) = 0;
	private:
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

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"
#include "BlumeCapel.h"
#include "XXZAntiferro.h"

const std::string outdir = "../out/";
const int myid = 0;

#include "marqov.h"
int main()
{
	//RegistryDB registry("./cfgs");

	// remove old output and prepare new one
	std::string command = "rm -r " + outdir;
	system(command.c_str());
	makeDir(outdir);
	

	// ---- live view ----
//		/*

		const int L = 30;
		RegularLattice lattice(L, 3);
		std::string outfile = outdir+"temp.h5";
		Marqov<RegularLattice, XXZAntiferro<double,double> > marqov(lattice, 1/0.46, outfile);
		marqov.init_hot();
		const int ncluster = L;
		const int nsweeps  = L/2; 
		marqov.wrmploop(50, ncluster, nsweeps);
		marqov.gameloop_liveview();
//		*/


	/*
	// ---- 2D/3D Ising testing section ----

	std::vector<int> nL = {8,12,16,24,32,48};

	const int dim = 3;

	const int    nbeta     = 12;
	const double betastart = 0.32;
	const double betaend   = 0.44;

	double betastep = (betaend - betastart) / double(nbeta);

	std::ofstream os;
	os.open("simplelog.dat");
	for (int i=0; i<nbeta; i++) os << betastart + i*betastep << endl;
	os.close();


	for (int j=0; j<nL.size(); j++)
	{
		const int L = nL[j];
		cout << endl << "L = " << L << endl << endl;
		makeDir(outdir+std::to_string(L));

		for (int i=0; i<nbeta; i++)
		{
			double currentbeta = betastart + i*betastep; 
			cout << "beta = " << currentbeta << endl;

			RegularLattice lattice(L, dim);
    
			std::string outfile = outdir+std::to_string(L)+"/"+std::to_string(i)+".h5";

//			Marqov<RegularLattice, Ising<int> > marqov(lattice, currentbeta, outfile);
			Marqov<RegularLattice, BlumeCapel<int> > marqov(lattice, currentbeta, outfile);

			marqov.init_hot();

			const int ncluster = L;
			const int nsweeps  = L/2; 

			marqov.wrmploop(2000, ncluster, nsweeps);
			marqov.gameloop(8000, ncluster, nsweeps);
		}
	}
	*/

//	/*

	// ---- O(3) testing section ----

	std::vector<int> nL = {8,12,16,24,32,48};

	int    nbeta     = 12;
	double betastart = 0.59;
	double betaend   = 0.71;

	double betastep = (betaend - betastart) / double(nbeta);

	std::ofstream os;
	os.open("simplelog.dat");
	for (int i=0; i<nbeta; i++) os << betastart + i*betastep << endl;
	os.close();

	for (int j=0; j<nL.size(); j++)
	{

		const int L = nL[j];

		cout << endl << "L = " << L << endl << endl;
		makeDir(outdir+std::to_string(L));

		for (int i=0; i<nbeta; i++)
		{
			double currentbeta = betastart + i*betastep; 
			cout << "beta = " << currentbeta << endl;

			RegularLattice lattice(L, 3);
		
			std::string outfile = outdir+std::to_string(L)+"/"+std::to_string(i)+".h5";
{
			Marqov<RegularLattice, Heisenberg<double,double> > marqov(lattice, currentbeta, outfile);
}
{
			Marqov<RegularLattice, Phi4<double,double> > marqov(lattice, currentbeta, outfile);
}
			Marqov<RegularLattice, XXZAntiferro<double,double> > marqov(lattice, currentbeta, outfile);

			marqov.init_hot();

			const int ncluster = L;
			const int nsweeps  = L/2; 

			marqov.wrmploop(300, ncluster, nsweeps);
			marqov.gameloop(500, ncluster, nsweeps);

		}
	}

//	*/

}
