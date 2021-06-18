#include <iostream>
#include <array>

//include the MARQOV library
#include "libmarqov/libmarqov.h"

//include the RegularLattice
#include "lattice/regular_hypercubic.h"

//include certain hamiltonian building blocks from the utilities
#include "hamiltonian/util/termcollection.h"
//include some predefined observables, e.g. the magnetization and the energy
#include "hamiltonian/util/observables.h"
//include initializers for now
#include "hamiltonian/util/initializers.h"

class MySimpleIsing
{
	public:
        //Define the StateVectors that this model will use. The Ising Model has Integers +-1, hence we go with plain ints
		typedef std::array<int, 1> StateVector;
        constexpr static int SymD = 1; //FIXME: Are there cases where SymD is different from the length of the state vector?
        // We need to get an initializer that defines how to initialize our state vectors.
		template <typename RNG>
		using MetroInitializer = Ising_Initializer<StateVector, RNG>;

		//  ----  Parameters  ----
		double J;
		const std::string name;// Every Hamiltonian MUST have a name. Required for the HDF5 output

		//  ----  Hamiltonian terms  ---- 
        std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};

        
		MySimpleIsing(double J) : J(J), name("MySimpleIsing"), obs_e(*this){}
		~MySimpleIsing() {delete interactions[0];}

		//  ----  Observables ----
		Magnetization  obs_m;
		Energy<MySimpleIsing>  obs_e;

        decltype(std::make_tuple(obs_m, obs_e)) observables = {std::make_tuple(obs_m, obs_e)};
};

using namespace std;
using namespace MARQOV;

int main()
{
    //initialize a lattice
    RegularHypercubic mylatt(8, 2); //2D 8x8 lattice
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp(".");
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);
    
    // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
    double beta = 4;
    double J = -2.1;
    auto hp = make_tuple(beta, J);
    
    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    
    //execute the Core routines of MARQOV.
    auto mysim = makeCore<RegularHypercubic, MySimpleIsing>(args);
    mysim.init();
    mysim.gameloop();
}
