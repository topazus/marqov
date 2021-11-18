#include <iostream>
#include <array>
#include <tuple>

//include the MARQOV library
#include "libmarqov/libmarqov.h"

//include the RegularLattice
#include "lattice/regular_hypercubic.h"
//include the CSV graph
#include "lattice/graph_from_csv.h"

#include "hamiltonian/Heisenberg.h"
#include "hamiltonian/Ising.h"
#include "hamiltonian/Phi4.h"
#include "hamiltonian/AshkinTeller.h"
#include "hamiltonian/BlumeCapel.h"
#include "hamiltonian/BlumeEmeryGriffiths.h"
#include "hamiltonian/Potts.h"

using namespace std;
using namespace MARQOV;

/** @file libshowcase.cpp
 * This files defines the function calls that are accessible via SWIG in other languages.
 */

/** Create a single simulation
 * 
 * @tparam Lattice The type of the Lattice
 * @tparam Ham The type of the hamiltonian
 * @tparam HamParameter the tuple of Hamiltonian parameters
 * 
 * @param path the path where we write the data
 * @param mylatt a reference to a lattice
 * @param hp the tuple with the hamiltonian parameters
 * 
 * @return 0
 */
template <class Lattice, class Ham, typename HamParameter>
static int instantiatesim(const std::string& path, Lattice& mylatt, const HamParameter& hp)
{
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp(path);
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);

    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);

    //execute core
    auto mysim = makeCore<Lattice, Ham>(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
    return 0;
}

/** Create a single instance of the Ising Model
 * 
 * @see Ising.h
 * 
 * @param path The path where to store the data.
 * @param dim The dimensionality of the hypercubic lattice.
 * @param len The linear dimension.
 * @param beta The inverse temperature.
 * @param J the interaction energy.
 */
int pyIsing(std::string path, int dim, int len, double beta, double J)
{
	// Construct hypercubic regular lattice of linear length len and dimension dim
    RegularHypercubic mylatt(len, dim);

    // Set Hamiltonian parameters, the coupling J, and the inverse temperature beta
    auto hp = make_tuple(beta, J);

	// Create and run a single simulation
    return instantiatesim<RegularHypercubic, Ising<double>>(path, mylatt, hp);
}


int pyIsingGraph(std::string outpath, std::string graphpath, double beta, double J)
{
	// Load lattice from a graph specification
    GraphFromCSV graph(graphpath);

	// Set Hamiltonian parameters, the coupling J, and the inverse temperature beta
    auto hp = make_tuple(beta, J);

	// Create and run a single simulation
    return instantiatesim<GraphFromCSV, Ising<double>>(outpath, graph, hp);
}

/** Create a single instance of the Heisenberg Model
 * 
 * @see Heisenberg.h
 * 
 * @param path The path where to store the data.
 * @param dim The dimensionality of the hypercubic lattice.
 * @param len The linear dimension.
 * @param beta The inverse temperature.
 * @param J the interaction energy.
 */

int pyHeisenberg(std::string path, int dim, int len, double beta, double J)
{
    RegularHypercubic mylatt(len, dim);
    auto hp = make_tuple(beta, J);
    return instantiatesim<RegularHypercubic, Heisenberg<double, double> >(path, mylatt, hp);
}

int pyHeisenbergGraph(std::string outpath, std::string graphpath, double beta, double J)
{
    GraphFromCSV graph(graphpath);
    auto hp = make_tuple(beta, J);
    return instantiatesim<GraphFromCSV, Heisenberg<double, double>>(outpath, graph, hp);
}
 
int pyPhi4(std::string path, int dim, int len, double beta, double lambda, double mass)
{
     RegularHypercubic mylatt(len, dim);
     auto hp = make_tuple(beta, beta, lambda, mass);
 
     return instantiatesim<RegularHypercubic, Phi4<double, double> >(path, mylatt, hp);
}

int pyAshkinTeller(std::string path, int dim, int len, double beta, double J, double K)
{
     RegularHypercubic mylatt(len, dim);
     auto hp = make_tuple(beta, J, K);
     return instantiatesim<RegularHypercubic, AshkinTeller >(path, mylatt, hp);
}

int pyAshkinTellerGraph(std::string outpath, std::string graphpath, double beta, double J, double K)
{
     GraphFromCSV graph(graphpath);
     auto hp = make_tuple(beta, J, K);
     return instantiatesim<GraphFromCSV, AshkinTeller >(outpath, graph, hp);
}

int pyBlumeCapel(std::string path, int dim, int len, double beta, double J, double D)
{
     RegularHypercubic mylatt(len, dim);
     auto hp = make_tuple(beta, J, D);
     return instantiatesim<RegularHypercubic, BlumeCapel<int> >(path, mylatt, hp);
}

int pyBlumeCapelGraph(std::string outpath, std::string graphpath, double beta, double J, double D)
{
     GraphFromCSV graph(graphpath);
     auto hp = make_tuple(beta, J, D);
     return instantiatesim<GraphFromCSV, BlumeCapel<int> >(outpath, graph, hp);
}

int pyBEG(std::string path, int dim, int len, double beta, double J, double D, double K)
{
     RegularHypercubic mylatt(len, dim);
     auto hp = make_tuple(beta, J, D, K);
     return instantiatesim<RegularHypercubic, BlumeEmeryGriffiths<int> >(path, mylatt, hp);
}


template <class Lattice, class HamParms>
static int instantiatePotts(int q, const std::string& path, Lattice&& mylatt, const HamParms& hp)
{ 
     //MARQOV::Config stores a set of Monte Carlo related parameters
     MARQOV::Config mp(path);
     mp.setnmetro(5);
     mp.setncluster(8/2);
     mp.setwarmupsteps(500);
     mp.setgameloopsteps(3000);

     //execute core
     switch (q)
     {
         case 3:
             return instantiatesim<typename std::decay<Lattice>::type, Potts<3> >(path, std::forward<Lattice>(mylatt), hp);
             break;
         case 4:
             return instantiatesim<typename std::decay<Lattice>::type, Potts<4> >(path, std::forward<Lattice>(mylatt), hp);
             break;
         case 6:
             return instantiatesim<typename std::decay<Lattice>::type, Potts<6> >(path, std::forward<Lattice>(mylatt), hp);
             break;
         case 8:
             return instantiatesim<typename std::decay<Lattice>::type, Potts<8> >(path, std::forward<Lattice>(mylatt), hp);
             break;
         default:
             throw(std::string("[MARQOV::Showcase] Supported Potts(q) values are q=3,4,6, and 8. All others only from C++."));
             break;
     }
     return 0;
    
}

int pyPotts(int q, std::string path, int dim, int len, double beta, double J)
{
     RegularHypercubic mylatt(len, dim);
     auto hp = make_tuple(beta, J);
     return instantiatePotts(q, path, mylatt, hp);
}
