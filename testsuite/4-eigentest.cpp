#include <iostream>
#include <array>
#include <tuple>
#include <vector>

#include <eigen3/Eigen/Core>

auto dot(const Eigen::Matrix<double, 3, 1>& a, const Eigen::Matrix<double, 3, 1>& b)
{
    return a.dot(b);
}

//inject eigen and tuple_size compatibility into the standard namespace.
namespace std {
	template <
		typename _Scalar,
		int _Rows, int _Cols,
		int _Options,
		int _MaxRows, int _MaxCols
	> class tuple_size<
		Eigen::Matrix<
			_Scalar,
			_Rows, _Cols,
			_Options,
			_MaxRows, _MaxCols
		>
	> {
	public:
		static_assert(
			_Rows != Eigen::Dynamic and _Cols != Eigen::Dynamic,
			"tuple_size is only supported for fixed size matrices"
		);
		static const size_t value = _Rows * _Cols;
	};
};

//include the MARQOV library
#include "libmarqov.h"

//include the RegularLattice
#include "../src/lattice/regular_hypercubic.h"

//include some predefined observables, e.g. the magnetization and the energy
#include "../src/hamiltonian/util/observables.h"



// Define interaction term for the Heisenberg model
class MyHeisenberg_interaction
{
	public:
		MyHeisenberg_interaction(const double J) : J(J) {}
		Eigen::Vector3d get (const Eigen::Vector3d& phi) {return phi;};
        const double J;
};

class MyEigenHeisenberg
{
	public:
		
		// The spin dimension of the Heisenberg or O(3) model
		constexpr static int SymD = 3;

		// Define the state vector that this model will use.
		// For the Heisenberg model this is a three-dimensional unit vector
		typedef Eigen::Vector3d StateVector;

		// Parameters
		double J; // The coupling constant
		const std::string name; // every Hamiltonian MUST have a name, this is required for the HDF5 output

		// Hamiltonian terms
		// here this is only the canonical O(3) interaction, defined above
		std::array<MyHeisenberg_interaction*, 1> interactions = {new MyHeisenberg_interaction(J)};
        
		// Constructor
		MyEigenHeisenberg(double J) : J(J), name("MyEigenHeisenberg"), obs_e(*this){}
		~MyEigenHeisenberg() {delete interactions[0];}

		// Observables
		Magnetization obs_m;
		Energy<MyEigenHeisenberg>  obs_e;
		std::pair<Magnetization, Energy<MyEigenHeisenberg> > observables{obs_m, obs_e};
//        decltype(std::make_tuple(obs_m, obs_e)) observables = {std::make_tuple(obs_m, obs_e)};
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
	MARQOV::Config mp("out"); // output path 
	mp.setnmetro(5); // number of Metropolis sweeps per EMCS
	mp.setncluster(10); // number of Wolff updates per EMCS
	mp.setwarmupsteps(500); // number of EMCS for warmup
	mp.setgameloopsteps(3000); // number of EMCS for production

	// Set the Hamiltonian parameters, J, and the inverse temperature beta
    double beta = 0.66;
    double J = 1;
    auto hp = make_tuple(beta, J);

    // Let's create some parameters tuples, a temperature scan for the scheduler to work on

    // Prepare argument list, contains
	// 1) reference to the lattice
	// 2) Monte Carlo parameter object
	// 3) Hamiltonian parameters packed as tuple
    auto args = make_tuple(std::ref(mylatt), mp, hp);
    
    // Eexecute the Core routines of MARQOV.
    auto mysim = makeCore<RegularHypercubic, MyEigenHeisenberg>(args);
    mysim.init();
    mysim.wrmploop();
    mysim.gameloop();
}
