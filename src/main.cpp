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

const int myid = 0;

#include "marqov.h"
int main()
{
	//RegistryDB registry("./cfgs");



//	/*
	// ---- 2D Ising testing section ----

	int    nbeta     = 1;
	double betastart = 0.44;
	double betaend   = 0.57;

	double betastep = (betaend - betastart) / double(nbeta);

	std::ofstream os;
	os.open("simplelog.dat");
	for (int i=0; i<nbeta; i++) os << betastart + i*betastep << endl;
	os.close();

	for (int i=0; i<nbeta; i++)
	{
		double currentbeta = betastart + i*betastep; 
		cout << "beta = " << currentbeta << endl;

		RegularLattice lattice(12, 2);
    
		std::string outfile = std::to_string(i)+".h5";

		Marqov<RegularLattice, Ising<int> > marqov(lattice, currentbeta, outfile);

		marqov.init_hot();
//		marqov.init_cold();
marqov.visualize_state_2d();
marqov.debugloop(1);

marqov.visualize_state_2d();
//		marqov.warmuploop(50);
//		marqov.gameloop(150);
	}
//	*/

	/*

	// ---- O(3) testing section ----

	int    nbeta     = 30;
	double betastart = 0.55;
	double betaend   = 0.85;

	double betastep = (betaend - betastart) / double(nbeta);

	std::ofstream os;
	os.open("simplelog.dat");
	for (int i=0; i<nbeta; i++) os << betastart + i*betastep << endl;
	os.close();

	for (int i=0; i<nbeta; i++)
	{
		double currentbeta = betastart + i*betastep; 
		cout << "beta = " << currentbeta << endl;

		RegularLattice lattice(32, 3);
	
		std::string outfile = std::to_string(i)+".h5";

		Marqov<RegularLattice, Heisenberg<double,double> > marqov(lattice, currentbeta, outfile);

		marqov.init_hot();
		marqov.warmuploop(50);
		marqov.gameloop(100);
	}

	*/

}
