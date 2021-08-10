#include <iostream>
#include <array>
#include <tuple>

//include the MARQOV library
#include "libmarqov/libmarqov.h"

//include the RegularLattice
#include "lattice/regular_hypercubic.h"

#include "hamiltonian/Heisenberg.h"
#include "hamiltonian/Ising.h"
#include "hamiltonian/Phi4.h"

using namespace std;
using namespace MARQOV;

int pyising(std::string path, int len, double beta, double J)
{
    RegularHypercubic mylatt(len, 2); //2D len x len lattice
    
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp(path);
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);
    
    // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
    auto hp = make_tuple(beta, J);
    //Let's create some parameters, a temperature scan for the scheduler to work on.
    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);

    //execute core
    auto mysim = makeCore<RegularHypercubic, Ising<double> >(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
    return 0;
}

int pyHeisenberg(std::string path, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, 2); //2D len x len lattice
//MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);

     // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
     auto hp = make_tuple(beta, J);
     //Let's create some parameters, a temperature scan for the scheduler to work on.
     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, Heisenberg<double, double> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
    return 0;
}
 
int pyPhi4(std::string path, int len, double beta, double lambda, double mass)
{
     RegularHypercubic mylatt(len, 2); //2D len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
     auto hp = make_tuple(beta, beta, lambda, mass);
     //Let's create some parameters, a temperature scan for the scheduler to work on.
     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, Phi4<double, double> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
     return 0;
}
//}
