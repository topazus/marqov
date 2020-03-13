#include <array>
#include <vector>
#include <iostream>
#include <string>
#include "rndwrapper.h"
#include "regular_lattice.h"
#include "vectorhelpers.h"

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



template <class StateVector>
class OnSite
{
	public: 
	    const StateVector h;
	    virtual StateVector operator() (StateVector& phi) = 0;
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

template<class VecType>
inline typename VecType::value_type dot(const VecType& a, const VecType& b)
{
    typedef typename VecType::value_type FPType;
    return std::inner_product(begin(a), end(a), begin(b), 0.0);
}


// ---------------------------------------


	 const int myid = 0;

#include "Heisenberg.h"
#include "Ising.h"

template <class StateSpace>
class Observable {
    public:
        const std::string name;
        Observable(const std::string s) : name(s) {}
        virtual void measure(const StateSpace&) = 0;
};

template <class StateSpace, typename FPType = double>
class Magnetization : public Observable<StateSpace> {
    public:
        Magnetization() : Observable<StateSpace>("Magnetization") {}
        double measure(const StateSpace& statespace) {return 1.0;};
};


template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;

		// Constructor
		Marqov(Grid& lattice) : ham(),  grid(lattice), rng(0, 1), metro(rng) 
		{
		  	rng.seed(42);
		  	rng.seed(time(NULL));
		  	rng.set_integer_range(lattice.size());
		  	statespace = new typename Hamiltonian::StateVector[lattice.size()];
		}

		// Definition of an EMCS
	    void elementaryMCstep()
	    {
	        for(int i = 0; i < grid.size(); ++i)
	        {
                const int rsite = rng.i(); // choose random site -> Moved one level above
	            metropolisstep(rsite);
	        }
	    }
	    
	    void gameloop()
	    {
	        for (int i = 0; i < nstep; ++i)
	            elementaryMCstep();
	        
	        for(int j = 0; j < nobs; ++j)
	            obs[j]->measure(statespace);
					//improve me: consider that there might be reuse across observables!
	    }
	
	
		 void print_state()
		 {
             for(int dim = 0; dim < 3; ++dim) {
			for(int i = 0; i < grid.length; ++i)
			{
				for(int j = 0; j < grid.length; ++j)
				{
					double current = statespace[grid.length*i+j][dim];

					if (current > 0.1) cout << "+ ";
					else if (current < -0.1) cout << "- ";
					else cout << "0 ";
				}
				cout << endl;
			}
			std::cout<<std::endl;
             }
		 }
	
		 void init_cold()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = 1;
			}
		 }

		 void init_hot()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, double, 3>(rng);
			}
		 }
	
	private:


	// Single Metropolis update step statevectors on a lattice
	// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
	inline int metropolisstep(int rsite)
	{
	    StateVector& svold = statespace[rsite];
	    StateVector svnew = metro.newsv(svold);

	    // interaction part
	    double interactionenergydiff = 0;
	    for(int a = 0; a < ham.Nalpha; ++a)
	    {
	        auto nbrs = grid.getnbrs(a, rsite);
	        StateVector averagevector = {0};
	
	        for (int i = 0; i < nbrs.size(); ++i)
	        {
	            auto mynbr = nbrs[i];
	            auto myvec = ham.interactions[a](statespace[mynbr]);
	            averagevector = averagevector + myvec;
	        }
	        interactionenergydiff += ham.interactions[a].J * (dot(svnew - svold, averagevector));
	    }
	
	    // onsite energy part
	    double onsiteenergyold = 0;
	    double onsiteenergynew = 0;
	    for (int b = 0; b < ham.Nbeta; ++b)
	    {
	        onsiteenergyold += dot(ham.onsite[b]->h, ham.onsite[b]->operator()(svold));
	        onsiteenergynew += dot(ham.onsite[b]->h, ham.onsite[b]->operator()(svnew));
	    }

	    
	    // multi-site energy
	    double multisiteenergyold = 0;
	    double multisiteenergynew = 0;
	    for (int g = 0; g < ham.Ngamma; ++g)
	    {
	        multisiteenergynew += ham.multisite[g]->operator()(svnew, rsite, statespace);//FIXME: think about this...
	        multisiteenergyold += ham.multisite[g]->operator()(svold, rsite, statespace);//FIXME: think about this...
	        //forgot k_gamma
	    }
	    
	    double dE 	= interactionenergydiff 
	    			+ (onsiteenergynew - onsiteenergyold) 
				+ (multisiteenergynew - multisiteenergyold);
	
	
	    int retval = 0;
	    if ( dE <= 0 )
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    else if (rng.d() < exp(-ham.beta*dE))
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    
	    return retval;
	}


	void wolff()
	{
	}


	StateSpace statespace;
	Hamiltonian ham;
	Grid& grid;
	RND rng;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11

	// number of observables
	static constexpr uint nobs = 0;
	Observable<StateSpace>* obs[5];

	// number of EMCS
	static constexpr int nstep = 250;
};

int main()
{

	RegularLattice lattice(35, 2);
// 	Marqov<RegularLattice, Ising<int> > marqov(lattice);
    
    Marqov<RegularLattice, Heisenberg<double, double> > marqov(lattice);

	marqov.init_hot();
	marqov.print_state();
	marqov.gameloop();
	cout << endl << endl;
	marqov.print_state();
}
