#include <iostream>
#include <array>
#include <tuple>
#include <vector>

//include the MARQOV library
#include "libmarqov/libmarqov.h"

//include the RegularLattice
#include "lattice/regular_hypercubic.h"

//include some predefined observables, e.g. the magnetization and the energy
#include "hamiltonian/util/observables.h"


// Define interaction term for the Heisenberg model
class MyHeisenberg_interaction
{
	public:
		MyHeisenberg_interaction(const double J) : J(J) {}
		std::array<double, 3> get (const std::array<double, 3>& phi) {return phi;};
        const double J;
};

class MySimpleHeisenberg
{
	public:
		
		// The spin dimension of the Heisenberg or O(3) model
		constexpr static int SymD = 3;

		// Define the state vector that this model will use.
		// For the Heisenberg model this is a three-dimensional unit vector
		typedef std::array<double, SymD> StateVector;

		// Parameters
		double J; // The coupling constant
		const std::string name; // every Hamiltonian MUST have a name, this is required for the HDF5 output

		// Hamiltonian terms
		// here this is only the canonical O(3) interaction, defined above
		std::array<MyHeisenberg_interaction*, 1> interactions = {new MyHeisenberg_interaction(J)};
        
		// Constructor
		MySimpleHeisenberg(double J) : J(J), name("MySimpleHeisenberg"), obs_e(*this){}
		~MySimpleHeisenberg() {delete interactions[0];}

		// Observables
		Magnetization  obs_m;
		Energy<MySimpleHeisenberg>  obs_e;
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
	MARQOV::Config mp("./"); // output path 
	mp.setnmetro(5); // number of Metropolis sweeps per EMCS
	mp.setncluster(10); // number of Wolff updates per EMCS
	mp.setwarmupsteps(500); // number of EMCS for warmup
	mp.setgameloopsteps(3000); // number of EMCS for production

	// Set the Hamiltonian parameters, J, and the inverse temperature beta
    double beta = 0.66;
    double J = 1;
    auto hp = make_tuple(beta, J);

    // Let's create some parameters tuples, a temperature scan for the scheduler to work on

    // Prepare an array of parameter tuples
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    vector<decltype(args)> v;

	// Fill
    for(int j = 0; j < 13; ++j)
    {
        hp = make_tuple(beta+j*0.1, J);
        v.push_back(make_tuple(std::ref(mylatt), mp, hp));
    }

    // Set up the MARQOV schedular
    typedef typename GetSchedulerType<MySimpleHeisenberg, RegularHypercubic, decltype(args)>::MarqovScheduler SchedulerType;
    SchedulerType sched(100,1); //MARQOV makes use of all available threads by default.
    
	// Feed parameters to the scheduler which creates the simulations
    for(auto p : v)
		sched.createSimfromParameter(p);
	// Run
    sched.start();
}
