#include <iostream>
#include <array>
#include <tuple>
#include <vector>

// Our small helper library to read in strings from command line files.
#include "libmarqov/util/registry.h"

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
		MySimpleHeisenberg(double J) : J(J), name("Heisenberg"), obs_e(*this){}
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
    // We utilize a small helper library of ours, the registry which reads Windows .ini style files to read config files
    RegistryDB reg;
    try
    {
        reg.init("../src/config", "ini");
    }
    catch(Registry_Exception& re)
    {
        std::cout<<"[MARQOV::main] Configuration directory not found! Assuming you're starting this MARQOV binary for the first time!"<<std::endl;
        std::cout<<"WELCOME TO MARQOV!"<<std::endl;
        std::cout<<"[MARQOV::main] To get you going we will generate and populate a configuration directory locally under ./config"<<std::endl;
        makeDir("./config");
        const auto filename = std::string{"./config/select.ini"};
        if(!fileexists(filename))
        {
            std::ofstream select(filename);
            select<<"[General]"<<'\n'<<"Hamiltonian = Heisenberg"<<'\n'<<"[END]"<<std::endl;
            reg.init("./config", "ini");
        }
        else
        {
            std::cout<<"[MARQOV::main] "<<filename<<" already exists, but is not usable. I would overwrite its content, hence I'm terminating now"<<std::endl;
            throw;
        }
    }
    
    // With the registry available we can now start to read in parameters from config files. Note that these parameters get replicated in the final HDF5 containers.
    
    std::string name{"Heisenberg"};// We are a Heisenberg Model. Saves some typing.
    std::string fn{name + ".ini"};
    //let's get a couple of parameters from the registry
    int nreplicas  = reg.Get<int>(fn, name, "rep" );
    auto nclusteramp   = reg.Get<double>(fn, "MC", "nclusteramp");
    auto nclusterexp   = reg.Get<int>(fn, "MC", "nclusterexp");
    auto nmetro        = reg.Get<int>(fn, "MC", "nmetro");
    auto warmupsteps   = reg.Get<int>(fn, "MC", "warmupsteps");
    auto measuresteps  = reg.Get<int>(fn, "MC", "measuresteps");
    int L = reg.Get<int>(fn, name, "L" );// a single lattice size
    auto dim  = reg.Get<int>(fn, name, "dim" );
    
    //let's create a single lattice.
    RegularHypercubic latt(L, dim);
    
    // Set the Hamiltonian parameters, J, and the inverse temperature beta
    auto beta = reg.Get<std::vector<double> >(fn, "Heisenberg", "beta");//read in a set of inverse temperatures
    auto J    = reg.Get<std::vector<double> >(fn, "Heisenberg", "J");//read in a set of coupling strengths
    auto hp = cart_prod(beta, J);//Form the cartesian product of all betas and all Js, the hamiltonian parameters
    
    // prepare
    std::string outpath = ".out/"+std::to_string(L)+"/";// by default we dump into the current directory.
    
    // Set Monte Carlo parameters using MARQOV::Config
    MARQOV::Config mp(outpath);// output path 
    mp.setnmetro(nmetro);// number of Metropolis sweeps per EMCS
    mp.setncluster(int(nclusteramp*pow(L,nclusterexp)));// number of Wolff updates per EMCS
    mp.setwarmupsteps(warmupsteps);// number of EMCS for warmup
    mp.setgameloopsteps(measuresteps);// number of EMCS for production
    
    makeDir(mp.outpath);// A small utility function that helps us to create folders
    
    auto params = finalize_parameter(latt, mp, hp);//bundle the lattice, the marqov parameters and the hamiltonian parameters
    auto rparams = replicator(params, nreplicas);//duplicate everything for the amount of replicas
    
    // Instantiate the scheduler which waits for new threads.
    auto sched = makeScheduler<MySimpleHeisenberg, RegularHypercubic> (rparams[0]);// makeScheduler can figure out a lot from one set of parameters
    
    for(auto p : rparams) sched.createSimfromParameter(p);//submit parameter set to scheduler.
    sched.start();
}
