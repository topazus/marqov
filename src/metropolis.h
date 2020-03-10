#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <vector>
#include <cmath>
#include "rndwrapper.h"

// Single Metropolis update step statevectors on a lattice
// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
template<class MCState, class Hamiltonian>
inline int metropolisstep(MCState& mcstate, Hamiltonian& ham)
{
    typedef typename Hamiltonian::StateVector StateVector;
    StateVector s = {0};
    const int rsite = rn.i(); // choose random site -> Move one level above
    StateVector svold = mcstate[rsite];
    StateVector svnew = ham.creatersv(svold); //get a new randomized state vector
    
    double interactionenergydifference = 0;
    for(int a = 0; a < ham.Nalpha; ++a)
    {
        auto nbrs = ham.grid.getnbrs(a, rsite);
        StateVector averagevector();
        for (int i = 0; i < nbrs.size(); ++i)
        {
            averagevector += ham.interaction[a](nbrs[i]);
        }
        interactionenergydifference += ham.interaction[a].J * (dot(svnew - svold, averagevector));
    }
    
    double onsiteenergyold = 0;
    double onsiteenergynew = 0;
    for (int b = 0; b < ham.Nbeta; ++b)
    {
        onsiteenergyold += dot(ham.onsite[b].h, ham.onsite[b](svold));
        onsiteenergynew += dot(ham.onsite[b].h, ham.onsite[b](svnew));
    }
    
    double multisiteenergyold = 0;
    double multisiteenergynew = 0;
    for (int g = 0; g < ham.Ngamma; ++g)
    {
        multisiteenergynew += ham.multisite[g](svnew, rsite, mcstate);//FIXME: think about this...
        multisiteenergyold += ham.multisite[g](svold, rsite, mcstate);//FIXME: think about this...
        //forgot k_gamma
    }
    
    double energydiff = (onsiteenergynew - onsiteenergyold) + interactionenergydifference
    + (multisiteenergynew - multisiteenergyold);
    
    int retval = 0;
    if ( energydiff > 0 )
    {
        svold = svnew;
        retval = 1;
    }
    else if (rn.d() < exp(-ham.beta*energydiff))
    {
        svold = svnew;
        retval = 1;
    }
    
    return retval;
}
#endif