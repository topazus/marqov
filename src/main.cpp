// -------------------
// The MARQOV project
// -------------------------------------------------------------------
// A modern framework for classical spin models on general topologies
// -------------------------------------------------------------------

#include <iostream>
#include <array>
#include <tuple>
#include <vector>

// Our small helper library to read in strings from config files.
#include "libmarqov/util/registry.h"

// Include the MARQOV library
#include "libmarqov/libmarqov.h"

// Include startup helpers
#include "libmarqov/util/startup.h"

// Include the regular hypercubic lattice
#include "lattice/regular_hypercubic.h"

// Include some predefined observables, e.g. the magnetization and the energy
#include "hamiltonian/util/observables.h"



// Define the data type of the state vector
// For the Heisenberg model this is a three-dimensional unit vector
typedef std::array<double, 3> HeisenbergSpin;

// Define interaction term for the Heisenberg model
class HeisenbergInteraction
{
	public:
		// Constructor
		HeisenbergInteraction(const double J) : J(J) {}

		// Interaction strength
        const double J;

		// Interaction type
		std::array<double, 3> get (const HeisenbergSpin& phi) {return phi;};
};




// Define the Heisenberg Hamiltonian
class HeisenbergHamiltonian
{
	public:
		
		// The spin dimension of the Heisenberg or O(3) model
		constexpr static int SymD = 3; // must be set

		// Define the data type of the state vector
		// For the Heisenberg model this is a three-dimensional unit vector
		typedef std::array<double, SymD> StateVector; // must be set

		// Parameters
		double J; // The coupling constant
		const std::string name; // every Hamiltonian MUST have a name

		// Hamiltonian terms
		// here this is only the canonical O(3) interaction, defined above
		std::array<HeisenbergInteraction*, 1> interactions = {new HeisenbergInteraction(J)};
        
		// Constructor
		HeisenbergHamiltonian(double J) : J(J), name("Heisenberg"), obs_e(*this){}
		// Destructor
		~HeisenbergHamiltonian() {delete interactions[0];}

		// Observables
		Magnetization  obs_m;
		Energy<HeisenbergHamiltonian>  obs_e{*this};

		// Group observables
               std::pair<Magnetization, Energy<HeisenbergHamiltonian> > observables{obs_m, obs_e};
};


using namespace std;
using namespace MARQOV;

int main(int argc, char* argv[])
{
#ifdef MPIMARQOV
    int threadingsupport;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &threadingsupport);
	//FIXME: maybe we get by with one level less.
    if(threadingsupport < MPI_THREAD_SERIALIZED)
    {
        std::cout << "[MARQOV::main] Couldn't initialize MPI! ";
		std::cout << "Requested threading level not supported." << std::endl;
        return -1;
    }
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) 
	{
#endif
    print_welcome_message();
#ifdef MPIMARQOV
    }
#endif

	// We are a Heisenberg Model. Saves some typing...
    std::string name{"Heisenberg"};
    std::string fn{name + ".ini"};


	// Input
	// -----
    // We utilize a small helper library of ours, the registry which 
	// reads Windows .ini style files to read config files
    RegistryDB reg;
	check_registry_availability(reg, name);
	check_registry_file_exists(reg, name);
    // With the registry available we can now start to read in parameters from config files. 
	// Note that these parameters get replicated in the final HDF5 containers.
    

	// Geometry
	// --------
    // Instantiate a lattice
    int L = reg.Get<int>(fn, name, "L" ); // linear lattice size
    auto dim  = reg.Get<int>(fn, name, "dim" ); // dimension
    RegularHypercubic latt(L, dim);
    


	// Hamltonian parameters
	// ---------------------
    // Set the Hamiltonian parameters, the coupling J, and the inverse temperature beta
	// these are read from a file, using the registry module
    auto beta = reg.Get<std::vector<double>>(fn, "Heisenberg", "beta");
    auto J    = reg.Get<std::vector<double>>(fn, "Heisenberg", "J");

	// Form the cartesian product of all betas and all J's, hence spanning
	// our parameter space
    auto hp = cart_prod(beta, J);
    


	// Output directory
	// ----------------

    // Prepare output directory using helper script
	std::string outbasedir = "./out/";
    makeDir(outbasedir);
    std::string outpath = outbasedir+std::to_string(L)+"/";
    makeDir(outpath);
    


    // Monte Carlo parameters 
	// ----------------------

    // First, we import them from the registry
    int nreplicas  = reg.Get<int>(fn, name, "rep" );
    auto nclusteramp   = reg.Get<double>(fn, "MC", "nclusteramp");
    auto nclusterexp   = reg.Get<int>(fn, "MC", "nclusterexp");
    auto nmetro        = reg.Get<int>(fn, "MC", "nmetro");
    auto warmupsteps   = reg.Get<int>(fn, "MC", "warmupsteps");
    auto measuresteps  = reg.Get<int>(fn, "MC", "measuresteps");
	auto maxruntime    = reg.Get<double>(fn, "MC", "maxruntimehours")*std::chrono::hours{1};
	
	// Then, we store them into a MARQOV::Config object
    MARQOV::Config mp(outpath); // output path 
    mp.setnmetro(nmetro); // number of Metropolis sweeps per EMCS
    mp.setncluster(int(nclusteramp*pow(L,nclusterexp))); // number of Wolff updates per EMCS
    mp.setwarmupsteps(warmupsteps); // number of EMCS for warmup
    mp.setgameloopsteps(measuresteps); // number of EMCS for production
    


	// Schedule and Run
	// ----------------

	// Bundle the lattice, the MC parameters and the Hamiltonian parameters
    auto paramsets = finalize_parameter(latt, mp, hp);
    
	// Duplicate everything for the amount of replicas
    auto rparamsets = replicator(paramsets, nreplicas);
    
    // Instantiate the scheduler which waits for new threads.
	// (makeScheduler can figure out a lot from one set of parameters)
    auto sched = makeScheduler<HeisenbergHamiltonian, RegularHypercubic> (rparamsets[0]);
    sched.setmaxruntime(maxruntime);
	// Submit parameter sets to the scheduler
    for (auto p : rparamsets) sched.createSimfromParameter(p);

	// Run!
    sched.start();
#ifdef MPIMARQOV
    MPI_Finalize();
#endif
}
