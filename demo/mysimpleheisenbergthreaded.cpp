#include <iostream>
#include <array>
#include <tuple>
#include <vector>

//include the MARQOV library
#include "../src/libmarqov/libmarqov.h"

//include the RegularLattice
#include "../src/lattice/regular_hypercubic.h"

//include some predefined observables, e.g. the magnetization and the energy
#include "../src/hamiltonian/util/observables.h"


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

namespace MARQOV
{
    template <class Lattice>
    struct MARQOV::Wolff<MySimpleHeisenberg, Lattice>
    {
        std::vector<int> cstack = std::vector<int>(4096/sizeof(int), 0);///< the size of the stack is meant to be preserved across different cluster processes.
        
        template <class RNG, class StateSpace>
        inline int move(const MySimpleHeisenberg& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
        {
            typedef MySimpleHeisenberg Hamiltonian;
            typedef typename Hamiltonian::StateVector StateVector;
            constexpr static int SymD = Hamiltonian::SymD;
            
            int q = 0;
            auto& seed = statespace[rsite];
            cstack[q] = rsite;
            
            int rdir = rng.integer(SymD); // random Cartesian direction
            seed[rdir] *= -1;
            
            int clustersize = 1;
            int current = 0;
            
            // plain Heisenberg model has only one interaction term
            const auto gcpl = ham.interactions[0]->J; 
            
            while (q>=0)
            {
                current = cstack[q];
                q--;
                
                const auto nbrs = grid.nbrs(0, current);
                if(q + nbrs.size() > cstack.size()) 
                    cstack.resize(2*cstack.size());
                for (std::size_t i = 0; i < nbrs.size(); ++i)
                {
                    const auto currentnbr = nbrs[i];
                    StateVector& candidate = statespace[currentnbr];
                    
                    const auto lcpl = statespace[current][rdir] * candidate[rdir];
                    const auto cpl  = gcpl*lcpl;
                    
                    if (wolff_update_accepted(cpl, beta, rng))
                    {
                        q++;
                        clustersize++;
                        cstack[q] = currentnbr;
                        candidate[rdir] = -candidate[rdir];
                    }
                }
            }
            return clustersize;
        }
    };
}

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
    double J = -1;
    auto hp = make_tuple(beta, J);

    // Let's create some parameters tuples, a temperature scan for the scheduler to work on

    // Prepare an array of parameter tuples
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    vector<decltype(args)> v;

	// Fill
    for(int j = 0; j < 1; ++j)
    {
        hp = make_tuple(beta+j*0.1, J);
        v.push_back(make_tuple(std::ref(mylatt), mp, hp));
    }

    // Set up the MARQOV schedular
    auto sched = makeScheduler<MySimpleHeisenberg, RegularHypercubic>(args, 1);
	// Feed parameters to the scheduler which creates the simulations
    for(auto p : v)
		sched.createSimfromParameter(p);
	// Run
    sched.start();
system("rm  beta0.700000_0.h5");
}
