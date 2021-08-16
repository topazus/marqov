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
#include "hamiltonian/AshkinTeller.h"
#include "hamiltonian/BlumeCapel.h"
#include "hamiltonian/BlumeEmeryGriffiths.h"
#include "hamiltonian/Potts.h"

using namespace std;
using namespace MARQOV;

template <class Lattice, class Ham, typename HamParameter>
static int instantiatesim(const std::string& path, Lattice& mylatt, const HamParameter& hp)
{
//     Lattice mylatt(len, dim); //nD len x len lattice
    
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp(path);
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);
    
    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);

    //execute core
    auto mysim = makeCore<RegularHypercubic, Ham>(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
    return 0;
}

int pyIsing(std::string path, int dim, int len, double beta, double J)
{
    RegularHypercubic mylatt(len, dim); //nD len x len lattice
    // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
    auto hp = make_tuple(beta, J);
    return instantiatesim<RegularHypercubic, Ising<double>>(path, mylatt, hp);
}

int pyHeisenberg(std::string path, int dim, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     
    // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
    auto hp = make_tuple(beta, J);
    return instantiatesim<RegularHypercubic, Heisenberg<double, double> >(path, mylatt, hp);
}
 
int pyPhi4(std::string path, int dim, int len, double beta, double lambda, double mass)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
     auto hp = make_tuple(beta, beta, lambda, mass);
 
     return instantiatesim<RegularHypercubic, Phi4<double, double> >(path, mylatt, hp);
}

int pyAshkinTeller(std::string path, int dim, int len, double beta, double J, double K)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, K);
     return instantiatesim<RegularHypercubic, AshkinTeller >(path, mylatt, hp);
}

int pyBlumeCapel(std::string path, int dim, int len, double beta, double J, double D)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, D);
     return instantiatesim<RegularHypercubic, BlumeCapel<int> >(path, mylatt, hp);
}

int pyBEG(std::string path, int dim, int len, double beta, double J, double D, double K)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, D, K);
     return instantiatesim<RegularHypercubic, BlumeEmeryGriffiths<int> >(path, mylatt, hp);
}

int pyPotts(int q, std::string path, int dim, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J);
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);

     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     switch (q)
     {
         case 3:
             return instantiatesim<RegularHypercubic, Potts<3> >(path, mylatt, hp);
             break;
         case 4:
             return instantiatesim<RegularHypercubic, Potts<4> >(path, mylatt, hp);
             break;
         case 6:
             return instantiatesim<RegularHypercubic, Potts<6> >(path, mylatt, hp);
             break;
         case 8:
             return instantiatesim<RegularHypercubic, Potts<8> >(path, mylatt, hp);
             break;
         default:
             throw(std::string("[MARQOV::Showcase] Supported q values 3,4,6, and 8. All others only from C++."));
             break;
     }
     return 0;
}
