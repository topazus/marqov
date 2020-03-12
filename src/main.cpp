#include <array>
#include <vector>
#include <iostream>
#include "rndwrapper.h"
#include "regular_lattice.h"

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
std::array<double, 1> operator + (StateVector lhs,  StateVector rhs)
{
    std::array<double, 1> res;
    res[0] = lhs[0] + rhs[0];
    return res;
}

template <class StateVector>
std::array<double, 1> operator - (StateVector lhs,  StateVector rhs)
{
    std::array<double, 1> res;
    res[0] = lhs[0] - rhs[0];
    return res;
}

double dot(std::array<double, 1> lhs,  std::array<double, 1> rhs)
{
    return lhs[0] * rhs[0];
}

// ---------------------------------------




template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
	Ising_interaction()
	{
		this->J = -1;	// ferro
		// this->J = +1; 	// antiferro
	}
	StateVector operator() (StateVector& phi) {return phi;};
};



template <typename SpinType, typename MyFPType>
class Ising
{
public:
    constexpr static double beta = 1/2.26918;
    constexpr static int SymD = 1;
    typedef std::array<SpinType, SymD> StateVector;
    typedef MyFPType FPType;
    static constexpr uint Nalpha = 1;
    static constexpr uint Nbeta = 0;
    static constexpr uint Ngamma = 0;
    // requires pointers
    Interaction<StateVector>* interactions;
    OnSite<StateVector>* onsite[Nbeta];
    MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
    Ising()
    {
        interactions = new Ising_interaction<StateVector>();
    }

    StateVector creatersv(const StateVector& osv) 
	 {
        StateVector retval(osv);
        retval[0] = -retval[0];
        return  retval;
	}
};



	 const int myid = 0;
    RND rn(0, 1);



#include "metropolis.h"
template <class StateSpace>
class Observable {public: void measure(StateSpace&) {};};



template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
	    typedef typename Hamiltonian::StateVector* StateSpace;
	    Marqov(Grid& lattice) : ham(),  grid(lattice) {
	        statespace = new typename Hamiltonian::StateVector[lattice.size()];
	    }
	    void elementaryMCstep()
	    {
	        for(int i = 0; i < grid.size(); ++i)
	        {
	            metropolisstep<StateSpace, Hamiltonian> (statespace, ham, grid);
	        }
	    }
	    
	    void gameloop()
	    {
	        for (int i = 0; i < nstep; ++i)
	            elementaryMCstep();
	        
	        for(int j = 0; j < nobs; ++j)
	            obs[j].measure(statespace);
					//improve me: consider that there might be reuse across observables!
	    }
	
	
		 void print_state()
		 {
			for(int i = 0; i < grid.length; ++i)
			{
				for(int j = 0; j < grid.length; ++j)
				{
					int current = statespace[grid.length*i+j][0];
					if (current == 1) cout << "o ";
					else if (current == -1) cout << ". ";
					else cout << "  ";
				}
				cout << endl;
			}
		 }
	
		 void init_cold()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = 1;
			}
		 }
	
	private:
	    StateSpace statespace;
	    Hamiltonian ham;
	    Grid& grid;
	    static constexpr uint nobs = 0;
	    Observable<StateSpace> obs[5];
	    static constexpr int nstep = 250;
};

void wolff();


int main()
{
	rn.seed(42);
	rn.seed(time(NULL));

	RegularLattice lattice(40, 2);
	rn.set_integer_range(lattice.size());
	Marqov<RegularLattice, Ising<double, double> > marqov(lattice);

	marqov.init_cold();
	marqov.print_state();
	marqov.gameloop();
	cout << endl << endl;
	marqov.print_state();
}
