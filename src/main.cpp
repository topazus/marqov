#include <array>
#include <vector>
#include <iostream>
#include <string>
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


// ---------------------------------------

#include "Heisenberg.h"
#include "Ising.h"
#include "Phi4.h"

const int myid = 0;

#include "marqov.h"
#include <cstdlib>
int main()
{

	//RegistryDB registry("./cfgs");

	RegularLattice lattice(30, 2);
    
	std::string outfile = "output.h5";
	Marqov<RegularLattice, Ising<int> > marqov(lattice, 1/2.2, outfile);

	marqov.init_cold();
	cout << endl << "Equilibrating ... " << endl;
	marqov.gameloop(100);
	marqov.gameloop_liveview(500,1);

//	marqov.visualize_state_2d();
}
