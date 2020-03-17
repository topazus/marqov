#ifndef ISING_H
#define ISING_H
#include <array>
#include <tuple>
#include <string>
#include <functional>


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class IsingMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			double mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += statespace[i][0];
			}

			return mag/double(N);
		}
		IsingMag() : name("m") {}
};

// Squared magnetization
class IsingMagSquare
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			float mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += pow(statespace[i][0],2);
			}
			return mag/double(N);
		}
		IsingMagSquare() : name("m2") {}
};


// Fourth power of magnetization
class IsingMagFour
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			float mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += pow(statespace[i][0],2);
			}
			return mag/double(N);
		}
		IsingMagFour() : name("m4") {}
};



// ----------------------------------------------------------------------

template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
	Ising_interaction()
	{
		this->J = -1;	// +1 ferro, -1 antiferro
	}
	StateVector operator() (StateVector& phi) {return phi;};
};


template <class StateVector, class RNG>
class Ising_Initializer
{
	public:
		Ising_Initializer()   {}
		Ising_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			retval[0] = -retval[0];
			return retval;
		};
};





// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class Ising
{
	public:
		double beta;
		constexpr static int SymD = 1;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = Ising_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		Ising(double mybeta) : beta(mybeta) {	interactions[0] = new Ising_interaction<StateVector>(); }
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
	    
	
		// instantiate and choose observables
		IsingMag       obs_m;
		IsingMagSquare obs_m2;
		IsingMagFour   obs_m4;
		auto getobs()
		{
			return std::make_tuple(obs_m, obs_m2, obs_m4);
		}
	
};
#endif
