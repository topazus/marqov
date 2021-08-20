#include <iostream>
#include <array>

// include the MARQOV library
#include "../src/libmarqov/libmarqov.h"

// include the RegularLattice
#include "../src/lattice/regular_hypercubic.h"

// include certain hamiltonian building blocks from the utilities
#include "../src/hamiltonian/util/termcollection.h"
// include some predefined observables, e.g. the magnetization and the energy
#include "../src/hamiltonian/util/observables.h"

class MySimpleIsing
{
	public:
		// The spin dimension of the Ising model
        constexpr static int SymD = 1;

        // Define the state vector that this model will use. 
		// The Ising Model has integers spins +1/-1, hence we go with plain ints
		typedef std::array<int, SymD> StateVector;

		// Parameters
		double J; // The coupling constant
		const std::string name; // every Hamiltonian MUST have a name, this is required for the HDF5 output

		// Hamiltonian terms
		// here this is only a canonical two-body interaction
        std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};
        
		// Constructor
		MySimpleIsing(double J) : J(J), name("MySimpleIsing"), obs_e(*this){}
		~MySimpleIsing() {delete interactions[0];}

		// Observables
		Magnetization  obs_m;
		Energy<MySimpleIsing>  obs_e;
        decltype(std::make_tuple(obs_m, obs_e)) observables = {std::make_tuple(obs_m, obs_e)};
};

using namespace std;
using namespace MARQOV;

int main()
{
    // Initialize the lattice
	int L = 8;
	int dim = 2;
    RegularHypercubic mylatt(L, dim);

    // Set Monte Carlo parameters using MARQOV::Config
    MARQOV::Config mp("out"); // output filename: will be out.h5
    mp.setnmetro(5); // number of Metropolis sweeps per EMCS
    mp.setncluster(10); // number of Wolff updates per EMCS
    mp.setwarmupsteps(500); // number of EMCS for warmup
    mp.setgameloopsteps(3000); // number of EMCS for production
    
    // Set the Hamiltonian parameters, J, and the inverse temperature beta
    double beta = 0.440;
    double J = 1;
    auto hp = make_tuple(beta, J);
    
    // Prepare argument list, contains
	// 1) reference to the lattice
	// 2) Monte Carlo parameter object
	// 3) Hamiltonian parameters packed as tuple
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    
    // Eexecute the Core routines of MARQOV.
    auto mysim = makeCore<RegularHypercubic, MySimpleIsing>(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
}
