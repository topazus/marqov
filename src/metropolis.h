

// Single Metropolis update step statevectors on a lattice
// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
template<class MCState, Hamiltonian>
inline int metropolisstep(MCState& mcstate)
{
    StateVector s = {0};
    const int rsite = model.rn.i(); // choose random site -> Move one level above
    StateVector svnew = ham.creatersv(svold); //get a new randomized state vector
    double energydiff = 0;
    
    interactionenergy = 0
    for(int a = 0; a < ham.Nalpha; ++a)
    {
        vector<int> nbrs = grid.getnbrs(a, rsite);
        for (int i = 0; i < nbrs.size(); ++i)
        {
            averagevector += ham.interaction[a](nbrs[i]);
        }
        interactionenergydifference[a] += ham.interaction[a].J * (dot(svnew - svold, averagevector));
    }
    
    onsiteenergy = 0;
    for (int b = 0; b < ham.Nbeta; ++b)
    {
        onsiteenergyold += dot(ham.onsite[b].h, ham.onsite[b](svold));
        onsiteenergynew += dot(ham.onsite[b].h, ham.onsite[b](svnew));
    }
    
    multisiteenergy = 0;
    for (int g = 0; g < ham.Ngamma; ++g)
    {
        multisiteenergynew += ham.multisite[g](svnew, rsite, statespace);//FIXME: think about this...
        multisiteenergyold += ham.multisite[g](svold, rsite, statespace);//FIXME: think about this...
        //forgot k_gamma
    }
    
    energydiff = (onsiteenergynew - onsiteenergyold) + (interactionenergynew - interactionenergyold) + (multisiteenergynew - multisiteenergyold);
    
    retval = 0;
    if ( energydiff > 0 )
    {
        svold = svnew;
        retval = 1
    }
    else if (rn.d() < exp(-beta*energydiff))
    {
        svold = svnew;
        retval = 1
    }
    
    return retval;
}
