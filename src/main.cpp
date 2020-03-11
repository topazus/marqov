#include <array>
#include <vector>
#include <iostream>
#include "rndwrapper.h"
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

class RegularLattice
{
public:
//    friend class NArray_Iterator;
    typedef std::vector<int> value_type;
    RegularLattice(int l, int d) : length(l), dim(d), pows(dim) {
        numberatoms = 1;
        for(int i = 0; i < dim; ++i)
        {
            pows[i] = numberatoms;
            numberatoms *= length;
        }
    }
    value_type getnbrs(int a, int i) const {return this->operator[](i);}
    value_type operator[](int i) const
    {
        std::vector<int> temp(2*dim);
        //calculate neighbours for site i
        for(int j = 0; j < dim; ++j)
        {
            int pl = pows[j]*length;
            //positive additions
            int c = i + pows[j];
            //test positive additions for PBCs
            if (c >= (c/pl)*pl)
            {//PBC
                temp[2*j] = (i/pl)*pl + c % pl;
            }
            else
                temp[2*j] = c;
            //      std::cout<<i<<" "<<temp[2*j]<<std::endl;
            
            //negative additions
            c = i - pows[j];
            //test negative additions for PBCs
            if(c < (i/pl)*pl)
            {
                temp[2*j+1] = (i/pl)*pl + (c + pl) % pl;
            }
            else
                temp[2*j+1] = c;
            //      std::cout<<i<<" "<<temp[2*j+1]<<std::endl;
        }
        return temp;
    }
    std::size_t size() const {return numberatoms;}
    int length;
private:
    std::size_t numberatoms;
    int dim;
    std::vector<int> pows;
};




class Grid {
    std::vector<int> getnbr(int i);
};

/** This 
 */
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

template
<class StateSpace, class StateVector>
class MultiSite
{
public: 
    const double k;
    virtual double operator() (StateVector& sv, int svpos, StateSpace s) = 0;
private:
};

std::array<double, 1> operator + (std::array<double, 1> lhs,  std::array<double, 1> rhs)
{
    std::array<double, 1> res;
    res[0] = lhs[0] + rhs[0];
    return res;
}

std::array<double, 1> operator - (std::array<double, 1> lhs,  std::array<double, 1> rhs)
{
    std::array<double, 1> res;
    res[0] = lhs[0] - rhs[0];
    return res;
}

double dot(std::array<double, 1> lhs,  std::array<double, 1> rhs)
{
    return lhs[0] * rhs[0];
}


template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
    Ising_interaction()
    {
        this->J = -1;	// ferro
//		  this->J = +1; 	// antiferro
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
class Marqov {
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
//	 rn.seed(myid+time(NULL)+random_device{}());

    RegularLattice lattice(20, 2);
    rn.set_integer_range(lattice.size());
    Marqov<RegularLattice, Ising<double, double> > marqov(lattice);
	 marqov.init_cold();
	 marqov.print_state();
    marqov.gameloop();
	 cout << endl;
	 marqov.print_state();
}
