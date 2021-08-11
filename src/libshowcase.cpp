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
    auto mysim = makeCore<RegularHypercubic, Ising<int> >(args);
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
/*    
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp(path);
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);
    

    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);

    //execute core
    auto mysim = makeCore<RegularHypercubic, Ising<double> >(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
    return 0;*/
}

int pyHeisenberg(std::string path, int dim, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
//MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);

     // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
     auto hp = make_tuple(beta, J);
     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, Heisenberg<double, double> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
    return 0;
}
 
int pyPhi4(std::string path, int dim, int len, double beta, double lambda, double mass)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
     auto hp = make_tuple(beta, beta, lambda, mass);
     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, Phi4<double, double> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
     return 0;
}

int pyAshkinTeller(std::string path, int dim, int len, double beta, double J, double K)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, K);

     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, AshkinTeller>(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
     return 0;
}

int pyBlumeCapel(std::string path, int dim, int len, double beta, double J, double D)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, D);

     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic, BlumeCapel<int> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
     return 0;
}

int pyBEG(std::string path, int dim, int len, double beta, double J, double D, double K)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J, D, K);

     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     auto mysim = makeCore<RegularHypercubic,BlumeEmeryGriffiths<int> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
     return 0;
}

int pyPotts(int q, std::string path, int dim, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, dim); //nD len x len lattice
 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);
 
     // A section for setting our Hamiltonian parameters, and beta
     auto hp = make_tuple(beta, J);

     //prepare the arguments
     auto args = make_tuple(std::ref(mylatt), mp, hp);
 
     //execute core
     switch (q)
     {
         case 3:
         {
     auto mysim = makeCore<RegularHypercubic, Potts<3> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
         }
             break;
         case 4:
         {
     auto mysim = makeCore<RegularHypercubic, Potts<4> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
         }
             break;
         case 6:
         {
     auto mysim = makeCore<RegularHypercubic, Potts<6> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
         }
             break;
         case 8:
         {
     auto mysim = makeCore<RegularHypercubic, Potts<8> >(args);
     mysim.init();
     mysim.wrmploop();
     mysim.gameloop();
         }
             break;
     }
     return 0;
}
