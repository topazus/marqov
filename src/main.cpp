#include "metropolis.h"

class Grid {
    std::vector<int> getnbr(int i);
};

/** This 
 */
struct Interaction
{
    const double J;
    double operator(StateVector phi_i, StateVector phi_j)
    {
        return (D_x * phi_x, D_y * phi_y ); //diagonal, unequal coupling within the structure of a statevector
        
        return 1/abs(i -j);//could also be from a table, e.g Dij[i,j]
    }
private:
    double [] Dij;
    Grid& grid;
};

struct OnSite
{
    const double[] h;
    StateVector operator(StateVector phi);
private:
    double [] data;
};

struct MultiSite
{
    const double[] k;
    StateVector operator(StateSpace s);
private:
    double data;
};

class Hamiltonian {
    typedef something StateVector;
    const uint Nalpha;
    const uint Nbeta;
    const uint Ngamma;
    
    Interaction[Nalpha] interactions;
    OnSite[Nbeta] onsite;
    MultiSite[Ngamma] multisite;
    
    StateVector creatersv(StateVector old);
};

class Marqov {
    void elementaryMCstep()
    {
        for(int i = 0; i < grid.size(); ++i)
        {
            metropolisstep(statespace);
        }
    }
    
    void gameloop()
    {
        for (int i = 0; i < nstep; ++i)
            elementaryMCstep();
        
        for(int j = 0; j < nobs; ++j)
            obs[j].measure(statespace);//FIXME: consider that there might be reuse across observables!
    }
    
};

void wolff();

class Observable {};



int main()
{
    
}
