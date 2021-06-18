#include <iostream>
#include <array>
#include <tuple>
#include <vector>

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

//Let's define our own interaction for the Heisenberg model
class MyHeisenberg_interaction
{
	public:
		MyHeisenberg_interaction(const double myJ) : J(-myJ) {}
		std::array<double, 3> get (const std::array<double, 3>& phi) {return phi;};
        const double J;
};

class MySimpleHeisenberg
{
	public:
        //Define the StateVectors that this model will use. The Heisenberg Model has real three component vectors.
		typedef std::array<double, 3> StateVector;
        constexpr static int SymD = 3; //FIXME: Are there cases where SymD is different from the length of the state vector?
        // We need to get an initializer that defines how to initialize our state vectors.
		template <typename RNG>
		using MetroInitializer = NVector_Initializer<StateVector, RNG>;

		//  ----  Parameters  ----
		double J;
		const std::string name;// Every Hamiltonian MUST have a name. Required for the HDF5 output

		//  ----  Hamiltonian terms  ---- 
        std::array<MyHeisenberg_interaction*, 1> interactions = {new MyHeisenberg_interaction(J)};
        
		MySimpleHeisenberg(double J) : J(J), name("MySimpleHeisenberg"), obs_e(*this){}
		~MySimpleHeisenberg() {delete interactions[0];}

		//  ----  Observables ----
		Magnetization  obs_m;
		Energy<MySimpleHeisenberg>  obs_e;

        decltype(std::make_tuple(obs_m, obs_e)) observables = {std::make_tuple(obs_m, obs_e)};
};

using namespace std;
using namespace MARQOV;

int main()
{
    //initialize a lattice
    RegularHypercubic mylatt(8, 2); //2D 8x8 lattice
    
    //MARQOV::Config stores a set of Monte Carlo related parameters
    MARQOV::Config mp("./");
    mp.setnmetro(5);
    mp.setncluster(8/2);
    mp.setwarmupsteps(500);
    mp.setgameloopsteps(3000);
    
    // A section for setting our Hamiltonian parameters, J, and the inverse temperature beta
    double beta = 4;
    double J = -2.1;
    auto hp = make_tuple(beta, J);
    //Let's create some parameters, a temperature scan for the scheduler to work on.
    //prepare the arguments
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    
    vector<decltype(args)> v;
    for(int j = 0; j < 7; ++j)
    {
        hp = make_tuple(beta+j*0.01, J);
        v.push_back(make_tuple(std::ref(mylatt), mp, hp));
    }

    //set up the scheduler of MARQOV
    typedef typename GetSchedulerType<MySimpleHeisenberg, RegularHypercubic, decltype(args)>::MarqovScheduler SchedulerType;
    SchedulerType sched(1);//MARQOV makes use of all available threads by default.
    
    for(auto p : v) sched.createSimfromParameter(p, defaultfilter);//feed parameters to the scheduler which creates the simulations.
    sched.start();//GoGoGo!!
}
