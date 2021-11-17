// -------------------
// The MARQOV project
// -------------------------------------------------------------------
// A modern framework for classical spin models on general topologies
// -------------------------------------------------------------------

#include <iostream>
#include <array>
#include <tuple>
#include <vector>

// Our small helper library to read in strings from command line files.
#include "libmarqov/util/registry.h"

// Include the MARQOV library
#include "libmarqov/libmarqov.h"

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
		Energy<HeisenbergHamiltonian>  obs_e;

		// Group observables
        decltype(std::make_tuple(obs_m, obs_e)) observables = {std::make_tuple(obs_m, obs_e)};
};


void create_config_select(std::string HamiltonianName)
{
	const auto filename = std::string{"./config/select.ini"};
	std::ofstream select(filename);
	select << "[General]" << '\n' << "Hamiltonian = " << HamiltonianName << '\n' << "[END]" << std::endl;
	select.close();
}



using namespace std;
using namespace MARQOV;

int main()
{
	// We are a Heisenberg Model. Saves some typing...
    std::string name{"Heisenberg"};
    std::string fn{name + ".ini"};


    // We utilize a small helper library of ours, the registry which reads Windows .ini style files to read config files
    RegistryDB reg;
    try
    {
        reg.init("./config", "ini");
    }
    catch(Registry_Exception& re)
    {
        std::cout << std::endl << "WELCOME TO MARQOV!" << std::endl << std::endl;
        std::cout << "[MARQOV::main] No configuration directory found! ";
		std::cout << "Assuming you're starting this MARQOV binary for the first time!" << std::endl;
        std::cout << "[MARQOV::main] To get you going we will generate and populate ";
		std::cout << "a configuration directory locally under ./config" << std::endl;

        makeDir("./config");
        const auto filename = std::string{"./config/select.ini"};

        if(!fileexists(filename))
        {
			create_config_select(name);
            reg.init("./config", "ini");
        }
        else
        {
            std::cout << "[MARQOV::main] "<< filename <<" already exists, but is not usable.";
			std::cout << "I would overwrite its content, hence I'm terminating now"<<std::endl;
            throw;
        }

    }


	try
	{
    	int L = reg.Get<int>(fn, name, "L" );
	}
	catch(Registry_cfgfile_not_found_Exception& e)
	{
        std::cout<<"[MARQOV::main] Unable to find Heisenberg config! Generating new one in ./config/"<<fn<<std::endl;
        ofstream cfg("./config/" + fn);
        cfg << "[Heisenberg]\n" << "L = 12\n" << "rep = 4\n" << "dim = 3\n";
		cfg << "beta = 0.67,0.68\n" << "J = -1.0\n\n";
		cfg << "[MC]\n" << "nmetro = " << 2 << "\nnclusteramp = "<< 0.5 << "\nnclusterexp = " << 1;
		cfg << "\nwarmupsteps = " << 500 << "\nmeasuresteps = "<< 5000 << "\n";
		cfg << "[IO]\n" << "outdir = ./out\n" << "[END]" << std::endl;
        reg.init("./config");                                                                                                                                                                                                           
    } 
			
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
    
	// Submit parameter sets to the scheduler
    for (auto p : rparamsets) sched.createSimfromParameter(p);

	// Run!
    sched.start();
}
